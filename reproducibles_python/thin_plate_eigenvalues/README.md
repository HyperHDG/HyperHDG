# Eigenvalue problems on networks of thin plates

This directory contains the numerical experiments of the manuscript

> T. Chaumont-Frelet, M. Dwivedi, and A. Rupp,
> *Eigenvalue problems on geometric hypergraphs: singular limits and HDG approximation*.

Three plates of thicknesses `(3/2, 3/2, 3) eps` meet along a straight junction of cross section
`eps Y_N` that carries the mass density `rho = 4 mu / (39 eps)`. For `eps -> 0`, the eigenvalues of
the diffusion operator on the thin domain converge to those of a limit problem on the hypergraph
of three unit squares, whose junction carries the mass density `mu = 24 / (17 pi)`. The plates are
clamped at their far ends, and all other boundaries are natural. The experiments study

| study      | what                                                                     | rates (theory)          |
|------------|--------------------------------------------------------------------------|-------------------------|
| `hdg`      | HDG method for the limit problem, degrees 0 to 3, uniform refinement      | `h^(2p+1)`, `h^(p+1)`   |
| `thin`     | thin domains for `eps = 2^-3, ..., 2^-10` versus the limit problem       | `eps`, `eps^(1/2)`      |
| `combined` | thin-domain eigenvalues versus HDG eigenvalues (degree 1), all `eps`, `h` | `eps + h^(2p+1)`        |

for the clusters `lambda_1` (simple), `lambda_2 = lambda_3 = pi^2` (double), and
`lambda_4 = 17 pi^2 / 16` (simple). The eigenvalue rates and the eigenspace rates (gap in the
weighted mass space `H_0`) are listed in this order.

## Running the experiments

Requirements: a C++20 compiler, CMake 3.18 or newer, LAPACK, Python 3.9 or newer with NumPy and
SciPy, and the git submodules of HyperHDG (nanobind). Then,

```sh
git clone https://github.com/HyperHDG/HyperHDG.git
cd HyperHDG
git fetch origin <commit given in the manuscript>
git checkout <commit given in the manuscript>
git submodule update --init --recursive
python3 -m pip install numpy scipy cmake        # if not available otherwise
python3 reproducibles_python/thin_plate_eigenvalues/run_experiments.py
```

The C++ module is compiled on the first run (into `.hyperhdg-cache/`, which takes about a minute),
and all studies take about six minutes on a laptop. `--only hdg thin combined` selects studies, and
`--quick` runs the coarsest levels only (a few seconds, without writing tables). The results are
written to `output/thin_plate_eigenvalues/results/` (ignored by git), or to the directory given by
`--output`:

- `hdg.json`, `thin.json`, `combined.json`: all computed eigenvalues, errors, gaps, residuals, and
  the numbers of unknowns,
- `table_hdg_eigenvalue_error.tex`, `table_hdg_gap.tex`, `table_thin_eigenvalue_error.tex`,
  `table_thin_gap.tex`, `table_combined.tex`: the bodies of the tables in the manuscript,
- `environment.json`: the versions of Python, NumPy, and SciPy, the platform, and the run times.

The subdirectory `results/` of this directory contains the results of the manuscript (computed on
macOS with Apple clang), which are the reference. At the end of a run, the new results are compared
with the reference, file by file: identical files are reported as such, and for differing JSON
files the largest absolute deviation is printed. Repeated runs on the same machine give identical
results; other machines reproduce them up to rounding errors, which affect the eigenvalue errors of
the HDG method below about `1e-11` and the thin-domain eigenvalues at about `1e-6`. Since the
result files contain no run times or machine information, they can also be compared directly, e.g.
by `diff -r reproducibles_python/thin_plate_eigenvalues/results output/thin_plate_eigenvalues/results`
(which only reports `environment.json` and the files of `--quick` runs). To update the reference, run the script with
`--output reproducibles_python/thin_plate_eigenvalues/results`.

## Files

- `limit.py`: exact eigenpairs of the limit problem (by separation of variables).
- `plates.py`: the hypergraph of three plates (a `.geo` file) and the evaluation of the HDG
  eigenvalue errors and eigenspace gaps.
- `thin.py`: the thin domains. Since neither the cross section nor the coefficients depend on the
  direction of the junction, and both caps are natural boundaries, the eigenfunctions are
  `w(x, y) cos(k pi (z + 1/2))`, where `w` solves `-div(grad w) + k^2 pi^2 w = lambda rho w` on the
  T-shaped cross section. These problems are solved by the same HDG method (degree 5) on meshes
  that are graded towards the junction and towards the reentrant corners. The averaging map takes
  thickness averages on the plates (extended into the collar by reflection) and the average over
  the junction cross section. A second computation on refined meshes estimates the accuracy.
- `common.py`: compilation of the module, `.geo` output, assembly, eigensolver (Lanczos in shift and
  invert mode with diagonal equilibration, iterative refinement, a fixed start vector, and
  purification of the eigenvectors, since the mass matrix is singular), Gauss quadrature, and gaps.
- `run_experiments.py`: runs the studies, writes the results and tables, and compares them with
  the reference in `results/`.

The C++ part consists of `LocalSolver::DiffusionEigsWeighted`
(`include/HyperHDG/local_solver/diffusion_eigs_weighted_ldgh.hxx`), which provides the local
stiffness and mass matrices of the HDG method with eliminated flux (plate weights, densities, and
hypernode mass are read from the hyperedge properties of the domain file), and
`GlobalLoop::GeneralizedEigenvalue` (`include/HyperHDG/global_loop/generalized_eigenvalue.hxx`),
which assembles them and evaluates discrete functions at quadrature points.
