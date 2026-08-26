# Explicit python bindings (nanobind)

An alternative to the runtime cython compile machinery (`import/` + `cython/`): the user
writes a small module `.cxx` that spells out the template instantiation it wants, CMake
builds a plain extension module, python imports it. Both paths drive the very same C++
classes — they agree bit for bit, see *Correctness anchor* below.

Nothing here has to be built ahead of time: `python/hyperhdg.py` generates the CMake project
for such a module, builds it into a cache directory and imports the result, so an example is
just a python script.

## Quickstart

```sh
git submodule update --init --recursive        # nanobind lives in submodules/nanobind.git
python python/diffusion.py domains/cross.geo --degree 2 --refine 8
```

`python/diffusion.py` is the command-line counterpart of
`reproducibles_python/fiber_network_diffusion.py`: LDG-H diffusion on a hypergraph of
one-dimensional hyperedges, condensed system assembled as a COO triplet and solved with
scipy's CG. Graph, polynomial degree, refinement level and penalty come from the command
line — the (space dimension, degree) instantiation is compiled the first time it is asked
for (~3.5 s) and cached afterwards (~8 ms).

The whole C++ side of the example is this:

```cpp
using HDG = GlobalLoop::Elliptic<
  Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
  LocalSolver::Diffusion<1, 2, 4, TestParametersSinEllipt, double>, std::vector<double> >;

NB_MODULE(hdg_diffusion, m) { HyperHDG::bind_python<HDG>(m, "Diffusion_D3_P2"); }
```

`HyperHDG::bind_python<GlobalLoopT>(m, name)` (`include/HyperHDG/bind_python.hxx`) binds
whatever subset of the global-loop protocol the instantiation provides — every method is
`requires`-guarded, so an elliptic loop gets no `set_data`/`make_initial`, an eigenvalue loop
no `residual_flux`, and none of the loops needs binding code of its own. One header replaces
the `.pyx`/`.pxd` template pair per loop.

On the python side the class is reached by attribute when its name is written into the
source, and by subscript when the name is assembled at runtime:

```python
module = hyperhdg.load(code)
module.Diffusion_D3_P2                     # name known while writing the code
module[f"Diffusion_D{dim}_P{degree}"]      # name known only while running it
```

## Convergence

`python/convergence.sh` sweeps the example over degrees and refinement levels and prints the
observed rates — `diffusion_convergence_elliptic.py` from the shell, on a graph read from
file (the whole table below takes 12 s from a cold cache):

```
$ bash python/convergence.sh
python/diffusion.py on domains/cross.geo

 degree  refinement   unknowns          error    rate
      1           1          5    9.62210e-02       -
      1           2          9    2.63289e-02    1.87
      1           4         17    6.96677e-03    1.92
      1           8         33    1.79621e-03    1.96
      1          16         65    4.56260e-04    1.98
      1          32        129    1.14990e-04    1.99

      2           1          5    1.16024e-02       -
      2           2          9    1.62594e-03    2.84
      2           4         17    2.14504e-04    2.92
      2           8         33    2.75150e-05    2.96
      2          16         65    3.48301e-06    2.98
      2          32        129    4.38093e-07    2.99

      3           1          5    1.13812e-03       -
      3           2          9    7.85654e-05    3.86
      3           4         17    5.13481e-06    3.94
      3           8         33    3.27742e-07    3.97
      3          16         65    2.06932e-08    3.99
      3          32        129    1.29980e-09    3.99
```

i.e. h^(degree+1). Note that the manufactured solution of `TestParametersSinEllipt` is
u = sin(pi/2 x) with a right hand side taken from its second derivative along x, so it solves
the network problem only on hyperedges parallel to the x axis; on a graph with oblique edges
(`domains/simplex_1_2.geo`, a fiber network) the rate column plateaus — the discretisation
still converges, but to a different limit.

## Correctness anchor

`python/check_vs_cython.py` runs both bindings over a matrix of (domain, degree, refinement,
tau) and compares everything that crosses the language boundary — right hand side, condensed
system matrix, CG solution, error:

```
simplex_1_2.geo        P1 ref 1 tau 1.0  n      3  error 6.136910361895e-01  max deviation 0.000e+00
simplex_1_2.geo        P2 ref 1 tau 1.0  n      3  error 6.680898820286e-01  max deviation 0.000e+00
simplex_1_2.geo        P3 ref 2 tau 1.0  n      6  error 6.659025722979e-01  max deviation 0.000e+00
injection_test.geo     P2 ref 4 tau 2.0  n     21  error 1.217103173353e+00  max deviation 0.000e+00
simplex_1_3.geo        P1 ref 1 tau 1.0  n      4  error 5.242802967533e-01  max deviation 0.000e+00
simplex_1_3.geo        P3 ref 3 tau 0.5  n     16  error 5.818433126717e-01  max deviation 0.000e+00
cross.geo              P2 ref 8 tau 1.0  n     33  error 2.751495442787e-05  max deviation 0.000e+00

worst deviation over 7 configurations: 0.000e+00
```

## `hyperhdg.py`: compiling a module on demand

```python
import hyperhdg

module = hyperhdg.load('''
#include <HyperHDG/bind_python.hxx>
...
NB_MODULE(my_module, m) { HyperHDG::bind_python<HDG>(m, "Diffusion_D2_P4"); }
''')
hdg = module.Diffusion_D2_P4("domains/simplex_1_2.geo", 1.0)
```

`hyperhdg.py` is a plain module next to the scripts that use it, so an example in `python/`
just does `import hyperhdg`. From elsewhere, put the directory on the import path and point
the module at the source tree:

```sh
export PYTHONPATH=$PWD/python
export HYPERHDG_DIR=$PWD          # only needed away from the source tree
```

`load()` compiles and imports, `compile()` just returns the `.so` path. The string is the
*whole* module — real C++, not template-parameter fragments — and it is compiled by a
generated, human-readable CMake project. Everything lives in `./.hyperhdg-cache/<name>/`
(override with `cache_dir=`): the generated `module.cxx` and `CMakeLists.txt` next to the
cmake build directory `build/`, which also holds the resulting `.so`. A content fingerprint
over code and options skips the build when nothing changed; `force=True` rebuilds regardless.

The user stays in control of the build:

- `cmake_args=[...]` is passed straight to CMake (compiler, flags, generator, ...).
- HyperHDG is located via the `hyperhdg=` argument or the `HYPERHDG_DIR` env var — either a
  source tree (header include-dirs, nanobind from its submodule) or an install prefix
  (`find_package(HyperHDG CONFIG)`, linking `HyperHDG::HyperHDG`). Default: the source tree
  this file lives in.
- nanobind with an *installed* HyperHDG: the generated project uses
  `find_package(nanobind CONFIG)`; if the python `nanobind` package is importable in the
  running interpreter (`pip install nanobind`), its bundled cmake config is injected via
  `-Dnanobind_DIR` automatically, so no system-wide nanobind is needed.

Extension modules can never be re-initialized in a CPython process; `load()` works around
this the way the legacy machinery did: by default (`hash_name=True`) the `NB_MODULE` name is
suffixed with a content hash, so *edited code simply loads as a new module in the same
session* — variants coexist, switching back to earlier code returns the cached module
instantly, and all variants share one build tree.

A module built by CMake instead — an `.so` from `nanobind_add_module(...)`, imported the
ordinary way — is the same thing without the cache; nothing in `bind_python.hxx` depends on
how the module was produced.

## Design rationale: why this instead of the runtime cython compile

Both systems make the user spell C++ template parameters; the difference is where they land
and what checks them. The legacy path
(`reproducibles_python/fiber_network_diffusion.py`):

```python
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "File<1,2,std::vector,Point<2,double> >"
const.geometry        = "File<1,2,std::vector,Point<2,double> >"
const.node_descriptor = "File<1,2,std::vector,Point<2,double> >"
const.local_solver    = "Diffusion<1,2,4,TestParametersSinEllipt,double>"
const.cython_replacements = ["string", "string"]
const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
PyDP = HyperHDG.include(const)          # text-substitutes .pyx/.pxd templates,
hdg  = PyDP("domains/fiber_network_1000.geo")   # runs cython + a C++ compiler, imports the .so
```

Same amount of template spelling as the `using HDG = ...` above — but in compiler-checked,
clangd-completed C++ instead of unchecked python strings. Note the legacy config additionally
leaks binding internals into user code: `cython_replacements` are raw cython `.pxd` type
strings.

| Claim | Legacy runtime cython | Explicit nanobind |
|---|---|---|
| **User controls the build** | Flags frozen into `cmake_cython.cfg` at configure; `--std=gnu++20` and cython flags hardcoded (`import_cxx/compile_prep.py`); falls back to a hardcoded `g++-10` + `/usr/include/pythonX.Y` when the cfg is missing (`import_cxx/cmake.py`) | An ordinary CMake project: compiler, flags, sanitizers, generator — `cmake_args=[...]`, or a target in the main build |
| **Transparency** | Generated code compiled at import time into `build/shared_objects/mod<sha256>.so`; staleness decided by 7 timestamp rules (`import_cxx/compile_prep.py`) | The generated project is readable C++ and CMake in `.hyperhdg-cache/`; staleness is one content fingerprint, or ninja's depfiles for a committed module |
| **Where errors surface** | Raw compiler stderr at `import` time, mid-script/notebook | A python exception carrying the CMake log — or, for a committed module, once at build time with IDE support |
| **Maintenance surface** | ~718 LoC: 346 python machinery (`import/HyperHDG/import_cxx/`) + 372 LoC cython templates — one `.pyx`/`.pxd` pair *per global loop* (5 pairs: elliptic, parabolic, 3× eigenvalue) | One generic ~190-line header for all loops, ~210 lines of cache/CMake glue |
| **Extensibility** | A new capability of a loop is reachable from python only after its `.pxd` declarations and numbered `CyReplaceNN` substitution slots are extended | `requires`-guarded lambdas: a method a loop offers is bound, one it does not is skipped — new loop methods appear in python without touching binding infrastructure |
| **Data crossing** | All vectors by-value `std::vector` copies (see `cython/elliptic.pyx`) | `std::vector` casters today; nanobind `ndarray` allows zero-copy numpy↔`std::span` where it matters |
| **Environment coupling** | Asserts the *exact* python minor version baked in at configure (`import_cxx/cmake.py`); expects `build/cmake_cython.cfg`, a layout that predates CMake presets (preset builds write `build/<preset>/`), with a hardcoded CI-runner path as the only alternative | Any python ≥ 3.9 with `Development.Module`; ABI encoded in the `.so` name the normal way |

### Measured (this machine, 2026-08-26, gcc 16.2.1, python 3.14, `-O2`)

| Operation | Time |
|---|---|
| Legacy: cold `HyperHDG.include(const)` — one config, compile at import | 4.0 s |
| Legacy: warm re-run (cache hit, timestamp checks + import) | 0.01 s |
| New: cold `hyperhdg.load()` of one instantiation | 3.5 s |
| New: cache hit | 8 ms |

Comparable both cold and warm; what differs is what is on disk afterwards and who decides how
it was compiled.

### What the legacy path does genuinely well

Honest concessions:

- **On-demand instantiation from a notebook**: trying a new polynomial degree needs no C++
  file and no build step — `include()` compiles it right there. That UX is not given up here:
  `hyperhdg.load()` does the same, with a regular generated CMake project underneath instead
  of string-substituted cython templates.
- A config that was compiled once stays cached; for a fixed set of configs the day-to-day
  experience is fine.

### How the legacy pipeline holds up today

The design dates from ~2021 and served the project well: one `include()` call gives any of
five global loops without a build system in sight. Reading it today, most of its problems
are the natural cost of implementing a build system by hand inside a library — pieces that
CMake and nanobind now provide natively:

- **Two build trees that don't know about each other.** CMake freezes compiler and flags
  into `cmake_cython.cfg` in *its* binary dir; the python side reads
  `<repo>/build/cmake_cython.cfg` and writes artifacts to
  `<repo>/build/{cython_files,shared_objects}` — both hardcoded relative to the package
  location (`import_cxx/paths.py`). The two coincide only while everyone configures with
  `-B build`; with CMake presets (`build/<preset>/`) the cfg is no longer found and the
  machinery silently falls back to a hardcoded `g++-10`.
- **Inseparable from the source checkout.** `main_dir()` resolves templates, headers, cfg
  and artifacts relative to the package file; there is no `setup.py`; the CMake install
  rules ship only the C++ headers. The python interface can only be used in-tree. *New
  path:* modules are ordinary `.so`s importable from anywhere, and `hyperhdg.py` also
  targets an installed HyperHDG prefix via `find_package`.
- **Environment locks and hand-rolled staleness.** The exact python minor version is
  asserted twice (cfg and `cython_log.txt`), and rebuild decisions rest on seven timestamp
  rules approximating dependency tracking — including "is the `.so` older than
  `compile_prep.py` itself". *New path:* that bookkeeping is delegated to CMake/ninja,
  which do it precisely (depfiles) and were built for it.
- **Shell-out error handling.** Each compile step is `os.system(...)` + `assert`; failures
  surface as raw compiler stderr in the middle of a script or notebook run.
- **A template pair per global loop.** Every loop needs a hand-maintained `.pyx`/`.pxd`
  pair kept in sync with the C++ signatures, plus numbered `CyReplaceNN` slots.

None of this is a criticism of the original choice: hand-written per-class extension
modules, the 2021 alternative, would have been worse, and the string-config UX was ahead of
its time. The machinery simply predates the tools that now make its job unnecessary.

### Where this path can still improve

- **Coverage.** Only the elliptic diffusion example exists; parabolic and eigenvalue loops
  are already covered by the `requires`-guards in `bind_python.hxx` but have no example.
- **Data crossing still copies.** Vectors cross the boundary as by-value `std::vector` both
  ways; zero-copy `nb::ndarray` variants are possible but not implemented.
- **Installed-prefix mode is untested in anger.** Nobody installs HyperHDG today; the
  `find_package` branch of `hyperhdg.py` works by construction but hasn't run against a real
  install, and `find_dependency(nanobind)` in the exported HyperHDG config remains to be
  added.
- **No ctest.** The examples are run by hand; wiring them into CTest means deciding whether a
  test may invoke CMake itself.

## Files

- `include/HyperHDG/bind_python.hxx` — generic, `requires`-guarded binding helper.
- `python/hyperhdg.py` — compile/load a module source string (generated CMake project, cached).
- `python/diffusion.py` — the example: graph and parameters from the command line.
- `python/convergence.sh` — convergence sweep over degrees and refinement levels.
- `python/check_vs_cython.py` — bit-for-bit comparison against the legacy cython bindings.
