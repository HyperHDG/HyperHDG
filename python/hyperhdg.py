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
holds the resulting extension module (.so, or .pyd on Windows). With hash_name=True (the
default of load), every code variant has its own directory, so switching between variants
does not rebuild them. Every call runs `cmake --build`, so an edited HyperHDG header
recompiles the module like any other source change; a content fingerprint (code + generated
project + cmake args) skips the far slower configure step when nothing changed, and
force=True recompiles regardless.

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

import contextlib
import hashlib
import importlib.machinery
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

__all__ = ["compile", "load"]


class Module(ModuleType):
    """An extension module that can also be subscripted.

    Attribute access reads well when the class name is written into the source,
    `mod.Diffusion_D2_P1`; when the name is only known at runtime, subscripting says the
    same without a getattr: `mod[f"Diffusion_D{dim}_P{degree}"]`.
    """

    def __getitem__(self, name):
        try:
            return getattr(self, name)
        except AttributeError:
            raise KeyError(f"{name} is not defined by module {self.__name__}") from None

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
  {root}
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


def _write_if_changed(path, content):
    """Write only on a real change: an untouched mtime keeps the build system from
    recompiling code that did not move."""
    if path.is_file() and path.read_text() == content:
        return
    path.write_text(content)


@contextlib.contextmanager
def _locked(path):
    """Hold an exclusive lock on path, so that processes building the same module (e.g. parallel
    runs of a parameter sweep) wait for each other instead of writing the same build tree."""
    with open(path, "a+") as handle:
        if os.name == "nt":
            import msvcrt
            while True:
                try:  # LK_LOCK retries for about 10 seconds, then raises
                    handle.seek(0)
                    msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
                    break
                except OSError:
                    pass
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _find_extension(build, name):
    """The extension module built for name: its file name ends with a suffix the running
    interpreter can import (e.g. .cpython-313-darwin.so, .abi3.so, .pyd), and multi-config
    generators (Visual Studio, Xcode) put it into a subdirectory such as Release/."""
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        found = sorted(build.rglob(name + suffix))
        if found:
            return found[0]
    return None


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
    """Compile a module source string; return the path of the built extension module.

    Artifacts live in <cache_dir>/<module name>/ (default ./.hyperhdg-cache/<name>/); the
    module stays in its build/ subdirectory unless output_dir is given, then it is copied
    there. The build system is always asked -- it, not the fingerprint, knows about the
    HyperHDG headers module.cxx includes -- so only cmake's configure step is fingerprinted.
    force=True reconfigures and recompiles regardless.

    hash_name=True suffixes the NB_MODULE name with a content hash (the legacy cython
    machinery's trick): every code variant becomes a distinctly named extension module
    with its own build directory, so edited code stays importable within one interpreter
    session and switching between variants does not rebuild them. The default keeps the
    user's name for a stable artifact.

    The module is built for the running interpreter (sys.executable).
    """
    base = _module_name(code)
    mode, root = _find_hyperhdg(hyperhdg)
    cache = Path(cache_dir).resolve() if cache_dir else Path(".hyperhdg-cache").resolve()
    if mode == "installed":
        try:  # pip-installed nanobind ships its cmake config out of CMake's sight
            import nanobind
            cmake_args = [*cmake_args, f"-Dnanobind_DIR={nanobind.cmake_dir()}"]
        except ImportError:
            pass
    # Build for the running interpreter: CMake would otherwise pick the first Python it finds,
    # whose extension modules cannot be imported here (e.g. the system Python on macOS).
    cmake_args = ["-DCMAKE_BUILD_TYPE=Release", f"-DPython_EXECUTABLE={sys.executable}",
                  *cmake_args]

    fingerprint = hashlib.sha256(
        "\0".join([code, mode, str(root), *map(str, cmake_args)]).encode()).hexdigest()
    name = f"{base}_{fingerprint[:8]}" if hash_name else base
    if hash_name:
        code = re.sub(r"NB_MODULE\(\s*" + re.escape(base), f"NB_MODULE({name}", code, count=1)
    work = cache / name  # one work/build tree per (possibly hash-suffixed) module name
    build = work / "build"

    if mode == "source":
        cmakelists = _SOURCE_TREE_CMAKE.format(name=name, root=root.as_posix())
        if not (root / "submodules/nanobind.git/CMakeLists.txt").is_file():
            raise FileNotFoundError(f"nanobind submodule missing in {root} "
                                    "(git submodule update --init --recursive)")
    else:
        paths = f"PATHS {root.as_posix()}" if root else ""
        cmakelists = _INSTALLED_CMAKE.format(name=name, package_paths=paths)

    work.mkdir(parents=True, exist_ok=True)
    with _locked(work / "lock"):
        return _build(work, build, name, code, cmakelists, fingerprint, cmake_args, force,
                      verbose, output_dir)


def _build(work, build, name, code, cmakelists, fingerprint, cmake_args, force, verbose,
           output_dir):
    stamp = work / f"{name}.fingerprint"
    configured = (build / "CMakeCache.txt").is_file()
    reconfigure = force or not configured or not stamp.is_file() \
        or stamp.read_text() != fingerprint

    stamp.unlink(missing_ok=True)  # a failed build must not leave a valid stamp behind
    _write_if_changed(work / "module.cxx", code)
    _write_if_changed(work / "CMakeLists.txt", cmakelists)
    if force:  # a fresh mtime recompiles even unchanged code
        (work / "module.cxx").touch()
    if reconfigure:
        _run(["cmake", "-S", ".", "-B", "build", *cmake_args], work, verbose)
    # always: only the build system knows whether an edited HyperHDG header invalidates the
    # cached .so -- the fingerprint sees module.cxx and the cmake arguments, not the headers
    # it includes. A no-op ninja build costs milliseconds.
    _run(["cmake", "--build", "build", "--parallel"], work, verbose)

    so = _find_extension(build, name)
    if so is None:
        raise RuntimeError(f"build produced no extension module {name} in {build}")
    stamp.write_text(fingerprint)
    return _deliver(so, output_dir)


def load(code, *, force=False, hash_name=True, **kwargs):
    """Compile a module source string and import it.

    hash_name=True (default) makes edited code reloadable within one session: each variant
    imports under a content-hash-suffixed module name (extension modules can never be
    re-initialized in-process, so a fresh name per variant is the only way).

    The module is returned as a Module, i.e. mod.TheClass and mod["TheClass"] both work.
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
    module.__class__ = Module  # allow subscripting for names built at runtime
    sys.modules[name] = module
    return module
