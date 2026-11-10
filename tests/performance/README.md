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
      -DHYPERHDG_CHOLMOD=ON -DHYPERHDG_BUILD_DOC=OFF -DNOPYTHONTESTS=ON
cmake --build build-perf --target performance
cd build-perf && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 srun --cpu-bind=cores -n 1 \
  ctest -L performance -V
```

Notes:

- Release is required. Debug keeps `hy_assert` and is `-O2` at most.
- `ctest -L performance` runs every benchmark serially (`RUN_SERIAL`), with the working directory
  set to the source root. Each benchmark prints Markdown tables and writes nanobench's JSON to
  `output/performance/<name>.json`.
- `linear_solvers` needs `-DHYPERHDG_CHOLMOD=ON` (SuiteSparse ≥ 7, CMake config
  `SuiteSparse::CHOLMOD`), the option that provides `SparseLA::Cholmod` in
  `include/HyperHDG/sparse_la.hxx`. Without it, the benchmark is left out. If the cluster has no
  SuiteSparse module, use spack (`spack install suitesparse`) or conda-forge, and point
  `CMAKE_PREFIX_PATH` at it.
- `timoshenko_layouts` takes a domain file as its argument (run the executable directly, from the
  source root); without one it refines `domains/cross.geo` to 100 000 beams. A `.geo.h5` file
  needs HDF5 (on by default if found).
- nanobench warns about frequency scaling and turbo. You cannot turn those off on the cluster;
  exclusive nodes keep the spread small anyway. Watch the `err%` column.

## What they measure

| executable | question |
|---|---|
| `local_operators_d2`, `local_operators_d3` | Per hyperedge, by layout and local solver: what do prepare, apply, residual, streaming apply (no store) and assemble cost? Variants: Diffusion Condensed, Diffusion Factored, DiffusionTensor, and Diffusion with one operator per hyperedge (unshared, as with variable coefficients). Degrees 1–3. |
| `parameters` | What can a parameter struct take at runtime, and at what cost? Compiled functions, a runtime variable, a runtime function pointer, a runtime `std::function`, and runtime `dirichlet_nodes`. 2D degree 2, unshared operators. It also checks that all variants compute the same residual, and that changing the boundary at runtime changes the constraints. |
| `linear_solvers` | When does a sparse direct solver (CHOLMOD) pay off against unpreconditioned CG, matrix-free and on the assembled matrix? Problem-size ladders in 2D/3D, degrees 1–3. The break-even is the number of right-hand sides with the same operator. It also checks that the direct solution is that of CG. |
| `timoshenko_layouts` | Per beam of a large network, by layout of the Timoshenko beam solver: what do prepare, apply, residual and assemble cost, and how large is the store? BeamStiffness (the default) against Condensed and Factored, degrees 1–3. It also checks that the layouts agree. |

## Results on a laptop (indicative only)

AMD Ryzen 7 250, GCC 16.2, OpenBLAS 0.3.34 (serial), SuiteSparse 7.14.1, Release, one thread,
frequency scaling on. Repeat on the cluster before drawing conclusions about absolute numbers. The
ratios below are what these runs are for.

**Local operators** (ns per hyperedge, prepared store, shared operator unless "unshared"):

| | 2D p1 | 2D p3 | 3D p1 | 3D p2 | 3D p3 |
|---|---:|---:|---:|---:|---:|
| Condensed apply | 166 | 179 | 271 | 590 | 1 445 |
| Factored apply | 224 | 731 | 599 | 3 699 | 18 254 |
| DiffusionTensor apply | 239 | 629 | 353 | 1 359 | 2 122 |
| Condensed, unshared apply | 166 | 189 | 298 | 1 143 | 6 343 |
| Condensed apply, streaming | 1 050 | 19 576 | 12 350 | 285 461 | 2 906 568 |
| Factored apply, streaming | 689 | 9 321 | 4 450 | 96 007 | 762 353 |
| DiffusionTensor apply, streaming | 290 | 980 | 399 | 1 508 | 2 829 |
| Condensed residual | 250 | 631 | 663 | 7 694 | 37 999 |
| Condensed, unshared prepare | 2 266 | 28 313 | 19 264 | 293 286 | 2 140 179 |
| stored bytes per operator (Condensed / Tensor) | 1 280 / 320 | 8 192 / 960 | 10 752 / 456 | 69 984 / 856 | 270 336 / 1 400 |

("streaming": without `prepare`, an operator is built per hyperedge and sweep.)

Observations:

- In 2D, apply costs about 170 ns per hyperedge at every degree, so per-edge overhead dominates,
  not the local product. The overhead is building the hyperedge from the topology and geometry and
  computing its dof indices on every sweep. Caching the dof indices of each hyperedge (a gather
  list, as a GPU would need anyway) is the obvious next optimization. The same overhead makes
  matrix-free CG slower than CG on the assembled matrix for 2D p1 (see below).
- Condensed is the fastest prepared apply. Factored is 2–13× slower in 3D (growing with the
  degree) because it re-solves through the factors. DiffusionTensor's prepared apply is 1.5×
  slower than Condensed at 3D p3, but its store is 190× smaller.
- Without `prepare`, Condensed is the wrong layout: it builds S and K for every hyperedge in
  every sweep. Factored costs a quarter of that at 3D p3, and DiffusionTensor a thousandth.
- Unshared (one operator per hyperedge) costs memory bandwidth on apply: 4.4× at 3D p3, where each
  operator is 270 kB. It costs flops on prepare: about 2 ms per hyperedge at 3D p3, which is
  quadrature plus factorization.
- `residual` is far more expensive than `apply` (26× at 3D p3). `local_load` re-integrates the
  right-hand side and the Dirichlet data on every sweep. Caching the load in the store for
  time-independent data, or sum-factorizing the load quadrature, would remove most of it.
- `assemble` of an unshared store builds every block with `apply<n_loc>` on the identity (n_loc³).
  A Condensed operator could hand out S directly instead.

**Parameters** (2D p2, unshared, relative to compiled functions):

| | prepare | residual | streaming apply |
|---|---:|---:|---:|
| runtime variable (`static inline double`) | +8 % | +2 % | +9 % |
| runtime function pointer | +18 % | +39 % | +18 % |
| runtime `std::function` (capturing lambda) | +19 % | +43 % | +20 % |
| runtime `dirichlet_nodes` (`std::vector`) | +8 % | +3 % | +9 % |

What can be runtime:

- The local solvers never instantiate the parameter struct; they name its static members. Diffusion
  passes the functions as template arguments to the quadrature, so the *functions* must be static
  and fixed at compile time.
- What those functions compute may depend on runtime state: static data members, which can be
  numbers, function pointers, or `std::function`. A static `std::vector<unsigned int>` also works
  for `dirichlet_nodes` and `neumann_nodes`.
- Compile-time only: everything the solver branches on with `if constexpr` or `requires`, namely
  `constant_coefficients` (operator sharing, explicit elimination) and the array types.
- A runtime number costs at most what the compiler loses in constant folding: 8 % in prepare here.
  The runtime boundary shows the same 8 % although prepare never reads the boundary, so part of
  it is code layout, not the variable. An indirect call per quadrature point costs 20–40 % with a
  cheap coefficient (1 + x²). With a coefficient that costs anything (a `sin`, a table lookup),
  the call disappears in the noise. With shared operators, the coefficient is evaluated only once
  per distinct operator anyway.
- The catch is that static state is global. Two loops with the same parameter type share it, and
  changing it between `prepare` and `apply` silently mixes old and new data. Per-loop runtime data
  would belong in the solver's `Args`, which is forwarded to every call.

**CG vs CHOLMOD** (seconds; unpreconditioned CG to rtol 1e-10 from zero; all dofs, the constrained
ones with a diagonal entry; break-even = right-hand sides after which CHOLMOD's setup beyond
`prepare` is paid back):

| case | dofs | factor nnz | CG its | MF CG | assembled CG | assemble | analyze | factorize | solve | break-even |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2D p1 n64 | 16 640 | 0.59 M | 339 | 0.241 | 0.050 | 0.0036 | 0.0033 | 0.012 | 0.00063 | 1 |
| 2D p1 n256 | 263 168 | 13.3 M | 1 314 | 14.9 | 3.89 | 0.044 | 0.070 | 0.261 | 0.027 | 1 |
| 2D p3 n64 | 33 280 | 2.34 M | 596 | 0.471 | 0.295 | 0.0072 | 0.0085 | 0.032 | 0.0038 | 1 |
| 3D p1 n16 | 52 224 | 12.8 M | 108 | 0.165 | 0.125 | 0.021 | 0.114 | 0.258 | 0.011 | 2.5 |
| 3D p2 n8 | 15 552 | 3.81 M | 76 | 0.027 | 0.049 | 0.020 | 0.021 | 0.060 | 0.0027 | 4.2 |
| 3D p2 n12 | 50 544 | 20.4 M | 112 | 0.137 | 0.258 | 0.073 | 0.095 | 0.482 | 0.014 | 5.3 |

(The full ladder is printed by the executable; "assemble" includes handing the triplets to
CHOLMOD, which compresses them.)

- In 2D, CHOLMOD wins from the first right-hand side at every size tried: nested-dissection fill
  is O(n log n), and unpreconditioned CG needs O(1/h) iterations (1 314 at 263k dofs).
- In 3D, CG wins a single solve. The break-even grows with the mesh, from 1.6 to 5.3 right-hand
  sides on this ladder, because fill and factorization grow like n^{4/3} and n², while CG
  iterations grow only like n^{1/3}. Memory goes the same way: the 3D p2 n12 factor holds 20 M
  nonzeros (160 MB), against a 70 kB shared operator for matrix-free.
- So a direct solver pays off in 2D, and in 3D when many solves share one operator: time stepping
  with a fixed step, many load cases, or the inner solve of a preconditioner. It does not pay for
  nonlinear or time-dependent coefficients that change the operator every solve. A preconditioned
  CG (multigrid or the network preconditioners) moves every break-even up. That comparison is the
  next one to add.
- Matrix-free CG is *slower* than CG on the assembled matrix at low degree (2D p1: 3.8×), the
  per-edge overhead from above. It is faster from 3D p2 on (n12: 0.137 vs 0.258 s), where the
  local product dominates.

**Timoshenko beams** (ns per beam; the cross refined to 100 000 beams, unit stiffnesses):

| degree | layout | prepare | apply | residual | assemble | bytes per beam |
|---:|---|---:|---:|---:|---:|---:|
| 1 | Condensed | 5 145 | 217 | 368 | 529 | 3 456 |
| 1 | BeamStiffness | 758 | 76 | 124 | 323 | 768 |
| 1 | Factored | 2 538 | 941 | 957 | 2 864 | 9 360 |
| 2 | Condensed | 12 151 | 228 | 331 | 516 | 4 608 |
| 2 | BeamStiffness | 1 746 | 77 | 146 | 312 | 1 088 |
| 2 | Factored | 6 751 | 1 684 | 1 737 | 5 764 | 16 056 |
| 3 | Condensed | 20 332 | 207 | 434 | 504 | 5 760 |
| 3 | BeamStiffness | 3 087 | 76 | 162 | 314 | 1 408 |
| 3 | Factored | 11 390 | 2 382 | 2 382 | 8 746 | 24 480 |

- BeamStiffness' apply costs the same at every degree and is compute bound: 4.1 instructions per
  cycle, against 1.6 for Condensed, which waits for its 330–550 MB store. It reads 128 bytes per
  beam; the rest of an operator (the fluxes of loads) is read by the residual only.
- prepare is 7× faster than Condensed's (four small factorizations per beam instead of one of
  the whole local matrix), the residual 2–3×, the store 4–5× smaller. All layouts agree to
  round-off (1e-15, printed by the executable).
- On a mikado network of 109 561 beams with the stiffnesses of its fibers (the command is in
  `timoshenko_layouts.cxx`): apply 125 against 270 ns, prepare 0.8–3.1 against 5.3–20.7 µs,
  residual 190–240 against 420–530 ns per beam.
