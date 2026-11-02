"""What every reproducible shares: the precompiled modules, the linear solver with its fallback
chain, the result line and output file, the plot options, and the parameter loops with the
command line (--quick for a short run).

A script defines one test function over its parameters and hands it to main() with the ranges:

    import common
    hdg_module = common.module("diffusion_elliptic")          # reproducibles_python/modules/
    def diffusion_test(poly_degree, dimension, iteration):
      hdg = hdg_module[f"Diffusion_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
      ...
      common.report("diffusion_convergence_elliptic", error, poly_degree=poly_degree, ...)
    if __name__ == "__main__":
      common.main(diffusion_test, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))

The modules are built by CMake (reproducibles_python/CMakeLists.txt) into
<build>/reproducibles_python/; <build> is $HYPERHDG_BUILD_DIR or build/ next to the sources.
"""
import inspect
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))   # hyperhdg.py is a script next to the examples, not a package yet
import hyperhdg  # noqa: E402

# SciPy renamed the tolerance of the iterative solvers from tol to rtol in version 1.12.
_TOLERANCE = "rtol" if "rtol" in inspect.signature(spla.cg).parameters else "tol"

# column labels of the output files, as the cython-era scripts wrote them
LABELS = {"poly_degree": "Polynomial degree", "dimension": "Dimension", "iteration": "Iteration",
          "refinement": "Refinement", "theta": "Theta", "tau": "Tau", "edge_dim": "Edge dimension",
          "space_dim": "Space dimension", "approx_type": "Approximation type"}


def module(name):
  """The precompiled module reproducibles_python/modules/<name>.cxx, subscriptable by class name."""
  return hyperhdg.import_prebuilt(name, subdir="reproducibles_python")


def solve(A, rhs, tol=1e-13, solvers=("cg", "gmres", "bicgstab")):
  """Solve A x = rhs with the first solver; on failure fall back along the chain with the
  messages the cython-era scripts printed (the reproducibles check greps them)."""
  for k, name in enumerate(solvers):
    solution, info = getattr(spla, name)(A, rhs, **{_TOLERANCE: tol})
    if info == 0:
      return solution
    if k + 1 < len(solvers):
      print(f"{name.upper()} failed with a total number of {info} iterations. Trying {solvers[k + 1].upper()}!")
    else:
      print(f"{name.upper()} also failed with a total number of {info} iterations.")
  raise RuntimeError("Linear solvers did not converge!")


def report(name, error, **params):
  """Print the error of one run and append it to output/<name>.txt, one line per run."""
  print("Iteration: ", params.get("iteration", params.get("refinement")), " Error: ", error)
  os.makedirs("output", exist_ok=True)
  with open(f"output/{name}.txt", "a") as f:
    f.write("".join(f"{LABELS.get(key, key)} = {value}. " for key, value in params.items())
            + f"Error = {error}.\n")


def plot(hdg, solution, filename, time=None, scale="0.95"):
  """Plot a solution to output/<filename>.vtu with the options the reproducibles use."""
  hdg.plot_option("fileName", filename)
  hdg.plot_option("printFileNumber", "false")
  hdg.plot_option("scale", scale)
  if time is None:
    hdg.plot_solution(list(solution))
  else:
    hdg.plot_solution(list(solution), time)


def main(test, argv=None, **grids):
  """Run test(**params) over the product of the grids (keyword order = loop nesting, the last
  one innermost). --quick (or the legacy positional 'True') keeps only the first two values of the
  innermost grid, a smoke test of every class. A RuntimeError of one run is printed as ERROR and
  the loops go on, as before."""
  argv = sys.argv[1:] if argv is None else argv
  quick = "--quick" in argv or "True" in argv
  grids = {key: list(values) for key, values in grids.items()}
  if quick:
    last = next(reversed(grids))
    grids[last] = grids[last][:2]

  def loops(keys, params):
    if not keys:
      start = datetime.now()
      print("Starting time is", start)
      try:
        test(**params)
      except RuntimeError as error:
        print("ERROR: ", error)
      end = datetime.now()
      print("Program ended at", end, "after", end - start)
      return
    key = keys[0]
    for value in grids[key]:
      print(f"\n{LABELS.get(key, key)} is set to be {value}\n")
      loops(keys[1:], {**params, key: value})

  loops(list(grids), {})
