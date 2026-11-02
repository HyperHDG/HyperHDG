# Reproducibles

Convergence studies of HyperHDG's global loops and local solvers, one script per study, on
precompiled nanobind modules. Every script defines one test function over its parameters and
hands it to `common.main` with the parameter ranges; `common.py` holds everything they share.

```
cmake -S . -B build && cmake --build build --target reproducibles      # the modules, once
python3 reproducibles_python/diffusion_convergence_elliptic.py          # the full study
python3 reproducibles_python/diffusion_convergence_elliptic.py --quick  # two refinements only
ctest --test-dir build -R reproducible_                                 # all, quick, checked
```

## Layout

- `<study>.py`: the study. About thirty lines: the module, the test function (set up the
  problem, solve, report the error, plot), and `common.main(test, order=..., **ranges)`.
- `modules/<family>.cxx`: one nanobind module per family of instantiations (global loop, local
  solver, topology). It defines an alias template `(dim, degree) -> global loop` and binds a
  grid of them as `<Name>_D<dim>_P<degree>` through `modules/bind_grid.hxx`; the grid is one
  line, so extending a module to another dimension or degree is a matter of a number. CMake
  builds them (`reproducibles_python/CMakeLists.txt`, option `HYPERHDG_REPRODUCIBLES`, target
  `reproducibles`) into `<build>/reproducibles_python/`, where `common.module` finds them
  (`$HYPERHDG_BUILD_DIR`, else `build/` next to the sources).
- `common.py`: module lookup, the condensed operator and the linear solver chain (CG, GMRES,
  BiCGStab), the time loop of the parabolic studies, the eigenvalue helpers (the reduced
  operators for scipy's `eigsh` and Newton's method on the eigenpair), the plot options, the
  records, the convergence summary and the check.
- `parameters/*.hxx`: the test problems (coefficients, right-hand sides, exact solutions).
- `expected/quick/<study>.yaml`, `expected/full/<study>.yaml`: the records a quick and a full run
  are checked against; the full ones are the reference of the convergence tables.

## Records

Every run appends one record to `output/<study>.yaml` and prints it, parameters first,
results after:

```
- {poly_degree: 2, dimension: 2, iteration: 3, error: 8.49e-07}
```

The file is a YAML list; `yq '.[] | select(.poly_degree == 2)' output/<study>.yaml` or PyYAML
read it whole. At the end `main` prints, per series (all parameters but the refinement one), the
errors and the experimental orders of convergence `log2(e_k / e_k+1)`, next to the expected
order the study states.

## Checking

- `--quick` keeps only the first two refinements, so every class runs once on coarse meshes;
  the records are a subset of the full run's.
- `--check` compares the records with the expected ones of the mode (every result to
  `atol + rtol * |expected|`, tolerances set per study in `main`), and the last rate of every
  series with three or more refinements with the expected order (up to `common.ORDER_TOL`); a
  mismatch, a missing record or a run that raised makes the script exit 1.
- `--record` writes the expected records of the mode from the run. CI checks the quick ones on
  every push; when a change moves a number beyond its tolerance, the full sweeps are run by hand,
  the convergence tables checked, and both sets re-recorded in that change.

The ctests run `--quick --check`; `shell_scripts/reproducibles_test.sh` runs every study the
same way and greps the logs, which the weekly workflow does.
