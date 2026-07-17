"""On-demand compilation of explicit HyperHDG nanobind modules.

The convenience layer on top of the explicit binding path (see python/README.md): pass the
complete module source (a tiny .cxx that calls HyperHDG::bind_python<...>) as a string,
get a compiled extension module back. Unlike the legacy cython machinery this generates a
regular CMake project in a visible cache directory, so the user keeps full control (any
cmake_args, choice of HyperHDG tree/install) and everything that is compiled is on disk in
plain sight.

    import hyperhdg
    mod = hyperhdg.load('''
      #include <HyperHDG/bind_python.hxx>
      ...
      NB_MODULE(my_module, m) { HyperHDG::bind_python<...>(m, "..."); }
    ''')

Everything lives in <cache_dir>/<module name>/ (default ./.hyperhdg-cache/<name>/): the
generated module.cxx and CMakeLists.txt next to the cmake build directory build/, which also
holds the resulting .so. A content fingerprint (code + generated project + cmake args) skips
the build when nothing changed; force=True rebuilds regardless.

HyperHDG is located in this order:
  1. `hyperhdg=` argument -- path to a source tree (has include/HyperHDG) or an install
     prefix (has lib/cmake/HyperHDG, used via find_package).
  2. environment variable HYPERHDG_DIR (same two flavors).
  3. the source tree this file lives in (python/ -> repo root).
  4. plain find_package(HyperHDG CONFIG) from whatever CMake sees.

nanobind: a source tree brings its own (submodules/nanobind.git). For an installed HyperHDG
the generated project calls find_package(nanobind CONFIG); if the python nanobind package is
importable (pip/pacman), its bundled cmake config is injected via -Dnanobind_DIR
automatically, so `pip install nanobind` into the running interpreter is all it takes.
"""

import hashlib
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

__all__ = ["compile", "load"]

_SOURCE_TREE_CMAKE = """\
cmake_minimum_required(VERSION 3.18)
project({name} LANGUAGES CXX)

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(LAPACK REQUIRED)
add_subdirectory({root}/submodules/nanobind.git nanobind)

nanobind_add_module({name} module.cxx)
target_compile_features({name} PRIVATE cxx_std_20)
target_include_directories(
  {name}
  PRIVATE
  {root}/include
  {root}/experiments
  {root}/submodules
  {root}/submodules/tensor_product_polynomials.git/include
  {root}/submodules/tensor_product_chain_complex.git/include
)
target_link_libraries({name} PRIVATE ${{LAPACK_LIBRARIES}})
"""

_INSTALLED_CMAKE = """\
cmake_minimum_required(VERSION 3.18)
project({name} LANGUAGES CXX)

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind CONFIG REQUIRED)
find_package(HyperHDG CONFIG REQUIRED {package_paths})

nanobind_add_module({name} module.cxx)
target_compile_features({name} PRIVATE cxx_std_20)
target_link_libraries({name} PRIVATE HyperHDG::HyperHDG)
"""


def _module_name(code):
    match = re.search(r"NB_MODULE\(\s*(\w+)", code)
    if not match:
        raise ValueError("code contains no NB_MODULE(<name>, m) definition")
    return match.group(1)


def _find_hyperhdg(hyperhdg):
    """Resolve the HyperHDG location to ('source'|'installed', path-or-None)."""
    candidates = [(hyperhdg, True), (os.environ.get("HYPERHDG_DIR"), True),
                  (Path(__file__).resolve().parent.parent, False)]
    for cand, explicit in candidates:
        if cand is None:
            continue
        cand = Path(cand).resolve()
        if (cand / "include/HyperHDG").is_dir() and (cand / "submodules").is_dir():
            return "source", cand
        if (cand / "lib/cmake/HyperHDG").is_dir():
            return "installed", cand
        if explicit:
            raise FileNotFoundError(f"{cand} is neither a HyperHDG source tree nor an "
                                    "install prefix")
    return "installed", None  # hope find_package knows better


def _run(cmd, cwd, verbose):
    result = subprocess.run(cmd, cwd=cwd, text=True,
                            capture_output=not verbose)
    if result.returncode != 0:
        output = "" if verbose else f"\n{result.stdout}\n{result.stderr}"
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed{output}")


def _deliver(so, output_dir):
    if output_dir is None:
        return so
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / so.name
    shutil.copy2(so, target)
    return target


def compile(code, *, hyperhdg=None, cache_dir=None, output_dir=None, cmake_args=(),
            force=False, hash_name=False, verbose=False):
    """Compile a module source string; return the path of the built .so.

    Artifacts live in <cache_dir>/<module name>/ (default ./.hyperhdg-cache/<name>/); the
    .so stays in its build/ subdirectory unless output_dir is given, then it is copied
    there. force=True recompiles regardless of the cached fingerprint.

    hash_name=True suffixes the NB_MODULE name with a content hash (the legacy cython
    machinery's trick): every code variant becomes a distinctly named extension module, so
    edited code stays importable within one interpreter session. The default keeps the
    user's name for a stable artifact.
    """
    base = _module_name(code)
    mode, root = _find_hyperhdg(hyperhdg)
    cache = Path(cache_dir).resolve() if cache_dir else Path(".hyperhdg-cache").resolve()
    work = cache / base  # one work/build tree per base name; variants share it
    build = work / "build"

    if mode == "installed":
        try:  # pip-installed nanobind ships its cmake config out of CMake's sight
            import nanobind
            cmake_args = [*cmake_args, f"-Dnanobind_DIR={nanobind.cmake_dir()}"]
        except ImportError:
            pass
    cmake_args = ["-DCMAKE_BUILD_TYPE=Release", *cmake_args]

    fingerprint = hashlib.sha256(
        "\0".join([code, mode, str(root), *map(str, cmake_args)]).encode()).hexdigest()
    name = f"{base}_{fingerprint[:8]}" if hash_name else base
    if hash_name:
        code = re.sub(r"NB_MODULE\(\s*" + re.escape(base), f"NB_MODULE({name}", code, count=1)

    if mode == "source":
        cmakelists = _SOURCE_TREE_CMAKE.format(name=name, root=root.as_posix())
        if not (root / "submodules/nanobind.git/CMakeLists.txt").is_file():
            raise FileNotFoundError(f"nanobind submodule missing in {root} "
                                    "(git submodule update --init --recursive)")
    else:
        paths = f"PATHS {root.as_posix()}" if root else ""
        cmakelists = _INSTALLED_CMAKE.format(name=name, package_paths=paths)

    # one stamp per (possibly hash-suffixed) name: alternating between two code variants of
    # the same base cache-hits both ways
    stamp = work / f"{name}.fingerprint"
    so = next(iter(sorted(build.glob(f"{name}.*.so"))), None)
    if not force and so and stamp.is_file() and stamp.read_text() == fingerprint:
        return _deliver(so, output_dir)

    work.mkdir(parents=True, exist_ok=True)
    stamp.unlink(missing_ok=True)  # a failed build must not leave a valid stamp behind
    # write_text refreshes module.cxx's mtime, so force=True recompiles even unchanged code
    (work / "module.cxx").write_text(code)
    (work / "CMakeLists.txt").write_text(cmakelists)
    _run(["cmake", "-S", ".", "-B", "build", *cmake_args], work, verbose)
    _run(["cmake", "--build", "build", "--parallel"], work, verbose)

    so = next(iter(sorted(build.glob(f"{name}.*.so"))), None)
    if so is None:
        raise RuntimeError(f"build produced no {name}.*.so in {build}")
    stamp.write_text(fingerprint)
    return _deliver(so, output_dir)


def load(code, *, force=False, hash_name=True, **kwargs):
    """Compile a module source string and import it.

    hash_name=True (default) makes edited code reloadable within one session: each variant
    imports under a content-hash-suffixed module name (extension modules can never be
    re-initialized in-process, so a fresh name per variant is the only way).
    """
    so = compile(code, force=force, hash_name=hash_name, **kwargs)
    name = so.name.split(".")[0]
    # extension modules cannot be re-initialized within a process: return the cached import
    if name in sys.modules:
        loaded = sys.modules[name]
        if force or getattr(loaded, "__file__", None) != str(so):
            raise ImportError(f"module '{name}' is already loaded from {loaded.__file__} "
                              "and extension modules cannot be reloaded in-process; use "
                              "hash_name=True (default), a different NB_MODULE name, or a "
                              "fresh interpreter")
        return loaded
    spec = importlib.util.spec_from_file_location(name, so)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[name] = module
    return module
