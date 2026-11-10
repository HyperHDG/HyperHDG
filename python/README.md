#Explicit python bindings(nanobind)

An alternative to the runtime cython compile machinery (`import/` + `cython/`): the user
writes a small module `.cxx` that spells out the template instantiation it wants, CMake
builds a plain extension module, python imports it. Both paths drive the very same C++
classes.

Nothing here has to be built ahead of time: `python/hyperhdg.py` generates the CMake project
for such a module, builds it into a cache directory and imports the result, so an example is
just a python script.

## Quickstart

```sh
git submodule update --init --recursive        # nanobind lives in submodules/nanobind.git
python3 examples/diffusion_elliptic_nb_jit.py --dim 2 --degree 2
```

The examples are `examples/diffusion_elliptic_nb_jit.py` (the module compiled on first use for
any dimension and degree), `examples/diffusion_elliptic_nb_module.cxx` with
`examples/diffusion_elliptic_nb_precompiled.py` (the same module built by CMake), all three run
as ctests.

Domains may also be `.geo.h5` files (the HDF5 counterpart of `.geo`, read and written by
`include/HyperHDG/read_domain_hdf5.hxx`, with a VTKHDF view for ParaView; `domains/cross.geo.h5`
is an example), if CMake finds HDF5: through `HDF5_ROOT` or `CMAKE_PREFIX_PATH` (a conda
environment: `$CONDA_PREFIX`), or `h5cc` on `PATH`.

### Plotting as VTKHDF

With HDF5 found, `plot_option("fileEnding", "vtkhdf")` makes `plot_solution(x, t)` write one HDF5
file `output/<fileName>.vtkhdf` in the VTKHDF format: the mesh on the first call, then one time
step per call (`t` is the step's time; `fileNumber` counts the steps), which ParaView opens as a
time series. `plot_option("valuesSelect", "mag:0,1,2")` or `"i,j"` restricts the stored
components of the solution, `"none"` drops them; a geometry's extra data becomes
`CellData/properties`. The `.vtu` writer stays the default.

### `.geo` to `.geo.h5`, and inspecting `.h5` files

```sh
python python/geo2h5.py domains/cross.geo                 # -> domains/cross.geo.h5
python python/geo2h5.py domains/cross.geo out.geo.h5 --no-vtkhdf
```

The dimensions are read from the header and the matching `read_domain` instantiation is compiled
on demand (a few seconds per combination, then cached). HDF5 ships `h5dump` and `h5ls` (packages
`hdf5` or `hdf5-tools`, or `hdf5` in a conda environment); none of these prints data unless asked:

```sh
h5dump -n  domains/cross.geo.h5               # list groups, datasets and links
h5dump -H  domains/cross.geo.h5               # shapes and types of everything, no data
h5dump -A  domains/cross.geo.h5               # attributes (space_dim, hyperedge_dim, ...)
h5dump -d /domain/edges domains/cross.geo.h5  # one dataset with its data
h5dump -p -d /domain/points domains/cross.geo.h5   # plus chunking and compression
h5ls -r    domains/cross.geo.h5               # one line per object, with shapes
```

From Python, h5py does the same:

```python
import h5py
with h5py.File("domains/cross.geo.h5") as f:
    f.visititems(lambda name, obj: print(name, getattr(obj, "shape", "")))
    points, edges = f["domain/points"][:], f["domain/edges"][:]
```

ParaView reads the `/VTKHDF` view as an unstructured grid. It picks readers by file extension and
does not know `.geo.h5`, so choose "VTKHDF Reader" when it asks, or symlink the file to a
`.vtkhdf` name.

The C++ side of an example is the instantiation and one line of binding, e.g.
`examples/diffusion_elliptic_nb_module.cxx`.

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

## `hyperhdg.py`: compiling a module on demand

```python
import hyperhdg

module = hyperhdg.load('''
#include <HyperHDG/bind_python.hxx>
...
NB_MODULE(my_module, m) {
  HyperHDG::bind_python<HDG>(m, "Diffusion_D2_P4"); }
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

`load()` compiles and imports, `compile()` just returns the path of the extension module
(`.so`, or `.pyd` on Windows). The string is the
*whole* module — real C++, not template-parameter fragments — and it is compiled by a
generated, human-readable CMake project. Everything lives in `./.hyperhdg-cache/<name>/`
(override with `cache_dir=`): the generated `module.cxx` and `CMakeLists.txt` next to the
cmake build directory `build/`, which also holds the resulting module. The module is built for
the running interpreter (`-DPython_EXECUTABLE=sys.executable`). A content fingerprint over code
and options skips the build when nothing changed; `force=True` rebuilds regardless. Processes
that build the same module at the same time wait for each other (a lock file in its directory).

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
instantly, and every variant has its own build tree (so switching does not rebuild).

A module built by CMake instead — an `.so` from `nanobind_add_module(...)`, imported the
ordinary way — is the same thing without the cache; nothing in `bind_python.hxx` depends on
how the module was produced.

## Design rationale: why this instead of the runtime cython compile

Both systems make the user spell C++ template parameters; the difference is where they land
and what checks them. The legacy path, as the cython-era reproducibles wrote it:

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
hdg  = PyDP("domains/cross.geo")   # runs cython + a C++ compiler, imports the .so
```

Same amount of template spelling as the `using` in `examples/diffusion_elliptic_nb_module.cxx`
— but in compiler-checked, clangd-completed C++ instead of unchecked python strings. Note the
legacy config additionally leaks binding internals into user code: `cython_replacements` are raw
cython `.pxd` type strings.

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
  `<repo>/build/{
  cython_files, shared_objects}` — both hardcoded relative to the package
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

## Files

- `include/HyperHDG/bind_python.hxx` — generic, `requires`-guarded binding helper.
- `python/hyperhdg.py` — compile/load a module source string (generated CMake project, cached).
- `examples/diffusion_elliptic_nb_*` — the examples, just in time and precompiled.

## Prototype: one global loop over local operators

`include/HyperHDG/global_loop/generic.hxx` (`GlobalLoop::Generic`) is a prototype of a single
global loop (design note `doc/local_operator_design.md`): a scatter–gather engine whose only
knowledge of the physics is an immutable per-hyperedge *local operator* that the local solver
returns from `local_operator(edge, Kind, args)`. The operator provides `apply<k>` (k traces at
once; assembly is `apply<n_loc>` of the identity, so the matrix is the matrix-free operator by
construction) and `residual`; the loop hands it the solver and the hyperedge in every call. The
loop owns the store (`prepare` / `release`: one operator per slot, as the solver numbers them
with `operator_slot`, else one per hyperedge; without `prepare` an operator is built per
hyperedge and sweep; each public sweep decides once and runs a branch-free body), the constraint
index set (Dirichlet rows: identity in `apply`, masked in `assemble`), `n_triplets` + `assemble`
into caller arrays, and `reduce` (with any combination, `+` by default);
`GlobalLoop::errors(loop, x)` is a reduction with the solver's `error_sq`. `prepare` and the
sweeps are templates on the kind of operator (`Stiffness` by default, `loop.apply<Mass>(x, y)`),
with one store per kind; a solver with other kinds lists them as `Kinds`. Nothing in the loop
allocates global vectors: every entry point takes vector-like inputs and outputs (`std::vector`,
`std::span`). Serial; parallel sweeps are a separate step.

`LocalSolver::Diffusion` (`include/HyperHDG/local_solver/diffusion_ldgh.hxx`) implements the
protocol next to its legacy interface. It sets up the blocks of the local problem, A x = B λ + l,
flux = C x + D λ (`local_blocks`), and its layout, a template parameter, decides what the
operator keeps (`include/HyperHDG/local_solver/local_operator.hxx`): by default
`LocalOperator::Condensed`, S = C A⁻¹ B + D and K = C A⁻¹ (apply is S λ, the residual S λ + K l),
or `LocalOperator::Factored`, the LU of A with B, C, D, or `LocalOperator::Streaming`, which
stores nothing and sets up the local problem in every call (one empty operator for all
hyperedges), or a layout of one's own. Diffusion numbers the operators when the parameters
declare `constant_coefficients` and the geometry gives the congruence class of a hyperedge
(`congruence_class`, as `Geometry::UnitCube`; `Geometry::File` does not yet, so there every
hyperedge has its own operator): congruent hyperedges then share one. Operators are
trivially copyable (`SmallMat` copy/move defaulted), the store is one contiguous array.

Bindings (`bind_python.hxx`, `requires`-guarded as before): NumPy in and out. The prototype's
`apply(x, y)`, `residual(x, r)` and `assemble(rows, cols, vals)` write into caller-allocated
arrays; Elliptic's vector methods gained NumPy overloads (lists still work), so both loops cross
the boundary at the same cost.

Tested by `tests/diffusion_elliptic.cxx` (convergence ladder, direct solve against matrix-free
CG); the driver, sweeps and roundoff comparison against Elliptic that the prototype started with
are in the history (`git log -- python/operator_loop.py`).
