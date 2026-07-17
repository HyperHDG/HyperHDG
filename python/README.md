# Explicit python bindings (nanobind)

This directory holds HyperHDG's explicit python binding path: a tiny, committed `.cxx` file
spells out the template instantiations that get compiled, CMake builds a plain extension
module, python imports it. It coexists with (and for the workflows here, replaces) the
runtime cython compile machinery in `import/` + `cython/`.

## Quickstart

```sh
git submodule update --init --recursive        # nanobind lives in submodules/nanobind.git
cmake --preset openblas
cmake --build --preset openblas --target timowave_py
python python/timowave.py -deg 1 -nx 4 -nt 32 -T 1 -domain domains/single1.geo
```

The module lands in `build/openblas/python/`; `timowave.py` adds that to `sys.path` itself
(`-build` selects the preset). Output keys mirror the PETSc driver
(`experiments/timowave.cxx`) so the same parsing pipelines work on both.

**Adding an instantiation** is one line in `python/timowave_py.cxx`:

```cpp
HyperHDG::bind_python<HDGTimoWave<3, 2, Complex>>(m, "TimoWave4_P3S2");
```

then rebuild the target (measured below: ~12 s) and `hy.TimoWave4_P3S2` exists. The generic
helper `HyperHDG::bind_python<GlobalLoopT>(m, name)` in `include/HyperHDG/bind_python.hxx`
binds whatever subset of the global-loop protocol the instantiation provides
(`requires`-guarded per method), so it serves any global loop — no per-loop wrapper files.

## Design rationale: why explicit modules instead of the runtime cython compile

Both systems make the user spell C++ template parameters; the difference is where they land
and what checks them. The legacy path
(`reproducibles_python/diffusion_hypergraph_convergence_elliptic.py`):

```python
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "Cubic<3,3>"
const.geometry        = "UnitCube<3,3,double>"
const.node_descriptor = "Cubic<3,3>"
const.local_solver    = "Diffusion<3,2,6,HG<3>::TestParametersQuadEllipt,double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
PyDP = HyperHDG.include(const)          # text-substitutes .pyx/.pxd templates,
hdg  = PyDP([4, 4, 4])                  # runs cython + C++ compiler, imports the .so
```

The explicit path (`python/timowave_py.cxx`):

```cpp
using HDGTimoWave = GlobalLoop::Hyperbolic<Topology::File<1, 3>, Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2 * poly_deg, TestTimoWave4, double, stages>,
  std::vector<ScalarT>>;
HyperHDG::bind_python<HDGTimoWave<1, 1, double>>(m, "TimoWave4_P1S1");
```
```python
import timowave_py as hy
hdg = hy.TimoWave4_P1S1("domains/single1.geo", [tau, dt, 0.0])
```

Same amount of template spelling — but in a compiler-checked, clangd-completed, git-versioned
file instead of unchecked python strings. Note the legacy config additionally leaks binding
internals into user code: `cython_replacements` are raw cython `.pxd` type strings.

### Claims and evidence

| Claim | Legacy runtime cython | Explicit nanobind |
|---|---|---|
| **User controls the build** | Flags frozen into `cmake_cython.cfg` at configure; `--std=gnu++20` and cython flags hardcoded (`import_cxx/compile_prep.py`); falls back to a hardcoded `g++-10` + `/usr/include/pythonX.Y` when the cfg is missing (`import_cxx/cmake.py`) | Ordinary CMake target: preset, compiler, flags, sanitizers, `compile_commands.json` — everything the rest of the build has |
| **Transparency / reproducibility** | Generated code compiled at import time into `build/shared_objects/mod<sha256>.so`; staleness decided by 7 timestamp rules (`import_cxx/compile_prep.py`) | What ran is a build artifact of a committed source file; "which instantiation produced this figure" is answered by git |
| **Where errors surface** | Raw compiler stderr at `import` time, on every user's machine, mid-script/notebook | At build time, once, with IDE support; imports never compile |
| **Maintenance surface** | ~718 LoC: 346 python machinery (`import/HyperHDG/import_cxx/`) + 372 LoC cython templates — one `.pyx`/`.pxd` pair *per global loop* (5 pairs, elliptic/parabolic/3× eigenvalue) | One generic ~175-line header for all loops; per-module cost is the tiny `.cxx` (~40 LoC) |
| **Extensibility** | The Gauss/timowave work (stage interface, `Gauss::StageTime`, complex scalars, span-based entries) never got cython templates — each would require new `.pxd` declarations plus more `CyReplaceNN` substitution slots | The same work was bound by `requires`-guarded lambdas; complex-scalar and stage entries came for free (`TimoWave4_P2S2`) |
| **Data crossing** | All vectors by-value `std::vector` copies (see `cython/elliptic.pyx`) | `std::vector` casters today; nanobind `ndarray` allows zero-copy numpy↔`std::span` where it matters |
| **Environment coupling** | Asserts the *exact* python minor version baked in at configure (`import_cxx/cmake.py`); expects `build/cmake_cython.cfg`, a layout that predates CMake presets (preset builds write `build/<preset>/`), with a hardcoded CI-runner path as the only alternative | Any python ≥ 3.9 with `Development.Module`; ABI encoded in the `.so` name the normal way |
| **PETSc coupling** | n/a (legacy loops are petsc-free) | Deliberately petsc-free: same headers compiled without `HYPERHDG_PETSC`; `ldd` on the `.so` shows no libpetsc |

Correctness anchor: `timowave.py` reproduces the PETSc driver
(`build/{openblas,complex}/experiments/timowave`) **digit-exactly** on wave4
(`-nx 4 -nt 32 -T 1`, single1.geo) for both the real single-stage (`P1S1`: `e_abs
1.24856e+01`, `e_rel 1.32982e+00`, `e_trace 4.45699e-01`) and the complex two-stage Gauss
instantiation (`P2S2`: `e_abs 2.47296e+00`, `e_rel 2.63390e-01`, `e_trace 1.97688e-02`).

### Measured (this machine, 2026-07-17, gcc 16.1.1, python 3.14, `-O3 -march=native`)

| Operation | Time |
|---|---|
| Legacy: cold `HyperHDG.include(const)` — one elliptic config, compile at import | 10.4 s |
| Legacy: warm re-run (cache hit, timestamp checks + import) | 0.07 s |
| New: clean `--target timowave_py` build (nanobind lib + 2 instantiations, parallel) | 10.4 s |
| New: incremental rebuild after adding one instantiation line | 12.5 s |
| New: `import timowave_py` | 61 ms |

The one-config legacy compile and the whole-module clean build cost the same wall time; the
difference is *when and where* it is paid: once per build tree at build time, versus at import
time on every machine and for every new config (and again whenever a timestamp rule fires).
The legacy warm path and the new import are equivalent (both instant).

*(Benchmark note: the legacy path currently needs `ln -s openblas/cmake_cython.cfg
build/cmake_cython.cfg` to run at all under CMake presets — its expected build layout predates
them.)*

### What the legacy path does genuinely well

Honest concessions:

- **On-demand instantiation from a notebook**: trying a new polynomial degree needs no C++
  file and no build step — `include()` compiles it right there. In the explicit system the
  equivalent is one line in `timowave_py.cxx` plus a ~12 s rebuild: cheap, but it is a step,
  and it lives outside the notebook.
- A config that was compiled once stays cached; for a fixed set of configs the day-to-day
  experience is fine.

That one-call UX exists *on top of* this system — see `hyperhdg.py` below: the same notebook
ergonomics, with a regular generated CMake project underneath instead of string-substituted
cython templates. The explicit path composes up to the convenient one; extracting build
control and transparency out of the string-substitution path does not work in reverse.

## On-demand compilation: `hyperhdg.py`

For exploratory work, `python/hyperhdg.py` compiles a module source string on demand:

```python
import hyperhdg

mod = hyperhdg.load('''
#include <HyperHDG/bind_python.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include <HyperHDG/local_solver/timowave.hxx>
#include "timowave4.hxx"

using HDG = GlobalLoop::Hyperbolic<
  Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoWave<1, 3, 1, 2, TestTimoWave4, double, 1>,
  std::vector<double>>;

NB_MODULE(jit_demo, m) { HyperHDG::bind_python<HDG>(m, "TimoWaveP1S1"); }
''')
hdg = mod.TimoWaveP1S1("domains/single1.geo", [tau, dt, 0.0])
```

`load()` compiles and imports; `compile()` just returns the `.so` path. The string is the
*whole* module — real C++, not template-parameter fragments — and it is compiled by a
generated, human-readable CMake project in `./.hyperhdg-build/<name>/` (override with
`build_dir=`); the `.so` lands in the current folder (override with `output_dir=`).
Content-hash caching skips the build when code and options are unchanged (measured: cold
~11 s, cache-hit load ~30 ms in a fresh process).

The user stays in control of the build:

- `cmake_args=[...]` is passed straight to CMake (compiler, flags, generator, ...).
- HyperHDG is located via the `hyperhdg=` argument or `HYPERHDG_DIR` env var — either a
  source tree (header include-dirs, nanobind from its submodule) or an install prefix
  (`find_package(HyperHDG CONFIG)`, linking `HyperHDG::HyperHDG`). Default: the source tree
  this file lives in.

One caveat inherent to CPython: extension modules cannot be re-initialized in a process, so
recompiling *different* code under an already-imported module name needs a new `NB_MODULE`
name or a fresh interpreter (`load()` raises a clear `ImportError` for this case).

### What becomes deletable

Once the loops in active use have explicit modules, `import/HyperHDG/import_cxx/` (346 LoC),
`cython/` (372 LoC in 5 template pairs + cfgs), and `cmake_cython.cfg.in` + its CMake glue
carry no remaining function and could be removed — including the cython, python-dev and
version-pinning requirements they impose on every user environment. Whether and when is a
separate decision; nothing here depends on it.

## Files

- `include/HyperHDG/bind_python.hxx` — generic, `requires`-guarded binding helper.
- `python/timowave_py.cxx` — the wave4 module: instantiation list, one line each.
- `python/CMakeLists.txt` — `nanobind_add_module`; petsc-free on purpose (links only
  tpp/tpcc/LAPACK, not the `HyperHDG` interface target).
- `python/timowave.py` — minimal scipy port of the PETSc driver's serial wave4 path.
- `python/hyperhdg.py` — on-demand compile/load of module source strings (generated CMake
  project, cached).
- `python/demo.sh` — 5-minute walkthrough of the above (run from the repo root).
