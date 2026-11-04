# Performance tests

Benchmarks with [nanobench](https://github.com/martinus/nanobench) (vendored in `extern/`, MIT).
They are not part of CI. CI runners are shared virtual machines, so their timings say nothing, and a
Release build of these executables alone takes longer than the whole test job.

## Where to run them: a bwUniCluster node

Run them on a **bwUniCluster compute node** (or another dedicated machine), not on a laptop and not
on a login node. You want one exclusive node, one thread and a fixed core:

```bash
salloc --partition=cpu --nodes=1 --exclusive --time=00:30:00   # or sbatch the lines below
module avail 2>&1 | grep -iE 'gnu|gcc|openblas|suitesparse|cmake'   # module names change
module load compiler/gnu numlib/openblas devel/cmake                # e.g.; pick from the list
cmake -S . -B build-perf -G Ninja -DCMAKE_BUILD_TYPE=Release -DHYPERHDG_PERFORMANCE_TESTS=ON \
      -DHYPERHDG_BUILD_DOC=OFF -DNOPYTHONTESTS=ON
cmake --build build-perf --target performance
cd build-perf && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 srun --cpu-bind=cores -n 1 \
  ctest -L performance -V
```

Notes:

- Release is required. Debug keeps `hy_assert` and is `-O2` at most.
- `ctest -L performance` runs every benchmark serially (`RUN_SERIAL`), with the working directory
  set to the source root. Each benchmark prints Markdown tables and writes nanobench's JSON to
  `output/performance/<name>.json`.
- `linear_solvers` needs CHOLMOD (SuiteSparse ≥ 7, CMake config `SuiteSparse::CHOLMOD`). Without
  it, the target is left out with a warning. If the cluster has no SuiteSparse module, use spack
  (`spack install suitesparse`) or conda-forge, and point `CMAKE_PREFIX_PATH` at it.
- nanobench warns about frequency scaling and turbo. You cannot turn those off on the cluster;
  exclusive nodes keep the spread small anyway. Watch the `err%` column.

## What they measure

| executable | question |
|---|---|
| `local_operators_d2`, `local_operators_d3` | Per hyperedge, by layout and local solver: what do prepare, apply, residual, streaming apply (no store) and assemble cost? Variants: Diffusion Condensed, Diffusion Factored, DiffusionTensor, and Diffusion with one operator per hyperedge (unshared, as with variable coefficients). Degrees 1–3. |
| `parameters` | What can a parameter struct take at runtime, and at what cost? Compiled functions, a runtime variable, a runtime function pointer, a runtime `std::function`, and runtime `dirichlet_nodes`. 2D degree 2, unshared operators. It also checks that all variants compute the same residual, and that changing the boundary at runtime changes the constraints. |
| `linear_solvers` | When does a sparse direct solver (CHOLMOD) pay off against unpreconditioned CG, matrix-free and on the assembled CSR? Problem-size ladders in 2D/3D, degrees 1–3. The break-even is the number of right-hand sides with the same operator. |

## Results on a laptop (indicative only)

AMD Ryzen 7 250, GCC 16.2, OpenBLAS 0.3.34 (serial), SuiteSparse 7.14.1, Release, one thread,
frequency scaling on. Repeat on the cluster before drawing conclusions about absolute numbers. The
ratios below are what these runs are for.

**Local operators** (ns per hyperedge, prepared store, shared operator unless "unshared"):

| | 2D p1 | 2D p3 | 3D p1 | 3D p2 | 3D p3 |
|---|---:|---:|---:|---:|---:|
| Condensed apply | 163 | 173 | 270 | 604 | 1 346 |
| Factored apply | 215 | 733 | 601 | 3 728 | 18 456 |
| DiffusionTensor apply | 192 | 653 | 383 | 1 399 | 2 109 |
| Condensed, unshared apply | 153 | 177 | 301 | 1 224 | 5 940 |
| Condensed apply, streaming | 582 | 9 376 | 3 763 | 92 165 | 761 428 |
| DiffusionTensor apply, streaming | 246 | 1 084 | 434 | 1 583 | 2 800 |
| Condensed residual | 237 | 627 | 673 | 7 714 | 38 279 |
| Condensed, unshared prepare | 2 184 | 28 875 | 18 900 | 292 288 | 2 086 304 |
| stored bytes per operator (Condensed / Tensor) | 1 280 / 320 | 8 192 / 960 | 10 752 / 456 | 69 984 / 856 | 270 336 / 1 400 |

Observations:

- In 2D, apply costs about 150 ns per hyperedge at every degree, so per-edge overhead dominates,
  not the local product. The overhead is building the hyperedge from the topology and geometry and
  computing its dof indices on every sweep. Caching the dof indices of each hyperedge (a gather
  list, as a GPU would need anyway) is the obvious next optimization. The same overhead makes
  matrix-free CG slower than CSR CG for 2D p1 (see below).
- Condensed is the fastest prepared apply. Factored is 2–14× slower in 3D (growing with the degree) because it re-solves
  through the factors. DiffusionTensor's prepared apply is 1.6× slower than Condensed at 3D p3, but
  its streaming apply is 270× faster than Condensed's, and its store is 190× smaller.
- Unshared (one operator per hyperedge) costs memory bandwidth on apply: 4.4× at 3D p3, where each
  operator is 270 kB. It costs flops on prepare: about 2 ms per hyperedge at 3D p3, which is
  quadrature plus factorization.
- `residual` is far more expensive than `apply` (28× at 3D p3). `local_load` re-integrates the
  right-hand side and the Dirichlet data on every sweep. Caching the load in the store for
  time-independent data, or sum-factorizing the load quadrature, would remove most of it.
- `assemble` of an unshared store builds every block with `apply<n_loc>` on the identity (n_loc³).
  A Condensed operator could hand out S directly instead.

**Parameters** (2D p2, unshared, relative to compiled functions):

| | prepare | residual | streaming apply |
|---|---:|---:|---:|
| runtime variable (`static inline double`) | +1.0 % | +0.4 % | +0.4 % |
| runtime function pointer | +10 % | +40 % | +19 % |
| runtime `std::function` (capturing lambda) | +17 % | +41 % | +32 % |
| runtime `dirichlet_nodes` (`std::vector`) | −0.5 % | +6 % | −0.7 % |

What can be runtime:

- The local solvers never instantiate the parameter struct; they name its static members. Diffusion
  passes the functions as template arguments to the quadrature, so the *functions* must be static
  and fixed at compile time.
- What those functions compute may depend on runtime state: static data members, which can be
  numbers, function pointers, or `std::function`. A static `std::vector<unsigned int>` also works
  for `dirichlet_nodes` and `neumann_nodes`.
- Compile-time only: everything the solver branches on with `if constexpr` or `requires`, namely
  `constant_coefficients` (operator sharing, explicit elimination) and the array types.
- The cost of runtime numbers is nil. An indirect call per quadrature point costs 10–40 % with a
  cheap coefficient (1 + x²). With a coefficient that costs anything (a `sin`, a table lookup),
  the call disappears in the noise. With shared operators, the coefficient is evaluated only once
  per distinct operator anyway.
- The catch is that static state is global. Two loops with the same parameter type share it, and
  changing it between `prepare` and `apply` silently mixes old and new data. Per-loop runtime data
  would belong in the solver's `Args`, which is forwarded to every call.

**CG vs CHOLMOD** (seconds; unpreconditioned CG to rtol 1e-10 from zero; break-even = right-hand
sides after which CHOLMOD's setup beyond `prepare` is paid back):

| case | dofs | factor nnz | CG its | MF CG | CSR CG | assemble | analyze | factorize | solve | break-even |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2D p1 n64 | 16 128 | 0.60 M | 339 | 0.223 | 0.040 | 0.0065 | 0.0031 | 0.012 | 0.00056 | 1 |
| 2D p1 n256 | 261 120 | 13.7 M | 1 314 | 13.9 | 3.55 | 0.086 | 0.076 | 0.283 | 0.028 | 1 |
| 2D p3 n64 | 32 256 | 2.40 M | 596 | 0.458 | 0.221 | 0.014 | 0.0086 | 0.033 | 0.0036 | 1 |
| 3D p1 n16 | 46 080 | 12.7 M | 108 | 0.129 | 0.099 | 0.032 | 0.108 | 0.259 | 0.011 | 3.4 |
| 3D p2 n8 | 12 096 | 3.76 M | 76 | 0.026 | 0.036 | 0.022 | 0.016 | 0.063 | 0.0024 | 4.3 |
| 3D p2 n12 | 42 768 | 19.3 M | 112 | 0.128 | 0.219 | 0.083 | 0.084 | 0.460 | 0.013 | 5.4 |

(The full ladder is printed by the executable.)

- In 2D, CHOLMOD wins from the first right-hand side at every size tried: nested-dissection fill
  is O(n log n), and unpreconditioned CG needs O(1/h) iterations (1 314 at 261k dofs).
- In 3D, CG wins a single solve. The break-even grows with the mesh, from 1.5 to 5.4 right-hand
  sides on this ladder, because fill and factorization grow like n^{4/3} and n², while CG
  iterations grow only like n^{1/3}. Memory goes the same way: the 3D p2 n12 factor holds 19 M
  nonzeros (150 MB), against a 70 kB shared operator for matrix-free.
- So a direct solver pays off in 2D, and in 3D when many solves share one operator: time stepping
  with a fixed step, many load cases, or the inner solve of a preconditioner. It does not pay for
  nonlinear or time-dependent coefficients that change the operator every solve. A preconditioned
  CG (multigrid or the network preconditioners) moves every break-even up. That comparison is the
  next one to add.
- Matrix-free CG is *slower* than CSR CG at low degree (2D p1: 3.9×), the per-edge overhead from
  above. It is faster from 3D p2 on (n12: 0.128 vs 0.219 s), where the local product dominates.
