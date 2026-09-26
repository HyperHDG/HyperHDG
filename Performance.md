HyperHDG evaluates the condensed HDG operator matrix-free by default: every call of `trace_to_flux`
solves the local problems of all hyperedges. This needs little memory and works for every local
solver, but it can be slow if the operator is applied many times. The options below speed up
typical computations. They do not change results (apart from rounding errors, where stated).


# Assembling the condensed matrix: `sparse_stiff_mat`

The global loops `Elliptic`, `Parabolic`, `ShiftedInverseEigenvalue`, and `MassApproxEigenvalue`
can assemble the condensed matrix in coordinate format:

```python
import scipy.sparse as sp

col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()   # Elliptic
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(time)   # other loops: time or shift
A = sp.csr_matrix((vals, (row_ind, col_ind)), shape=(n, n))   # sums up duplicate entries
```

The local solvers `Diffusion`, `DiffusionIP`, `DiffusionPostprocess`, `DiffusionAdvectionReaction`,
`DiffusionParab`, `Bilaplacian`, and `BilaplacianParab` factorize their local matrix only once per
hyperedge for this; the other local solvers solve the local problem for every local unit vector.
The resulting matrix is the same in both cases.

## Option `factorize_once`

```python
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(factorize_once=True)
```

computes the local matrix only once for every combination of boundary conditions of the faces
and reuses it for all other hyperedges. This is several times faster, but **it is only correct if
all hyperedges have the same local matrix**, i.e., for uniform meshes (e.g. `UnitCube` geometry on a
`Cubic` topology) with coefficients that are constant in space. It is **wrong** if, for example, the
diffusion tensor is space dependent or the hyperedges differ in size or shape. If a script runs with
`const.debug_mode = True`, this assumption is checked for every hyperedge and the program stops with
an error message if it is violated.


# Caching the operator: `cache_operator`

```python
HDG_wrapper.cache_operator(True)
A = LinearOperator((n, n), matvec=HDG_wrapper.trace_to_flux)   # as before
```

After `cache_operator(True)`, the first call of `trace_to_flux` assembles the condensed matrix and all
further calls evaluate a sparse matrix–vector product. This is useful for iterative solvers (CG,
GMRES, ...) and time stepping, where `trace_to_flux` is called very often. On an Apple M-series
computer, CG for a 3D diffusion problem with polynomial degree 2 became 25 times faster, and the
time stepping of `reproducibles_python/diffusion_convergence_parabolic.py` 6.8 times faster. The
results agree with the matrix-free ones up to rounding errors. The price is the memory for the matrix
(for a hyperedge of dimension d with n degrees of freedom per face, (2 d n)² values per hyperedge).

The arguments are `cache_operator(enable=True, factorize_once=False, time_independent=False)`
(`shift_independent` for `ShiftedInverseEigenvalue`). The matrix is assembled again whenever
`trace_to_flux` is called with another time (or shift), unless `time_independent=True` is set, and
after a refinement. `factorize_once` has the same meaning and restrictions as above.
`NonlinearEigenvalue` does not provide the cache, since its operator is not linear.


# Parallelization with OpenMP

The assembly of matrices (and the product with cached matrices) can be parallelized with OpenMP.
This is switched off by default and enabled by configuring HyperHDG with

    $ cmake .. -DHYPERHDG_OPENMP=ON

GCC and Clang on Linux as well as MinGW-w64 GCC on Windows support OpenMP out of the box. Apple
Clang on macOS needs the `libomp` library, e.g. from Homebrew:

    $ brew install libomp
    $ cmake .. -DHYPERHDG_OPENMP=ON -DOpenMP_ROOT=$(brew --prefix libomp)

The number of threads is set by the environment variable `OMP_NUM_THREADS`. The assembled matrices
do not depend on the number of threads. The matrix-free `trace_to_flux` is not parallelized, since
neighboring hyperedges write to the same entries; use `cache_operator` instead.
