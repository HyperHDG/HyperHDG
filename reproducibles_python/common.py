"""What every reproducible shares: the precompiled modules, the linear solver with its fallback
chain, the time loop, the eigenvalue helpers, the plot options, and the parameter loops with the
command line, the YAML records of the results, their convergence rates and their check against
the expected records.

A script defines one test function over its parameters and hands it to main() with the ranges:

    import common
    NAME = "diffusion_convergence_elliptic"
    hdg_module = common.module("diffusion_elliptic")          # reproducibles_python/modules/
    def diffusion_test(poly_degree, dimension, iteration):
      hdg = hdg_module[f"Diffusion_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
      ...
      common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
    if __name__ == "__main__":
      common.main(diffusion_test, order=lambda p: p["poly_degree"] + 1,
                  poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))

Every run appends one record to output/<NAME>.yaml and prints it,
    - {poly_degree: 1, dimension: 2, iteration: 3, error: 0.0012}
parameters first, results after (the file is a YAML list; read it whole with yq or PyYAML). At
the end main() prints, per series (all parameters but the refinement one), the errors and the
experimental orders of convergence log2(e_k / e_k+1).

Command line:  --quick    only the first two refinements (every class once, on coarse meshes);
                          the time-stepping scripts also take far fewer steps (common.QUICK)
               --check    compare the records with expected/{quick,full}/<NAME>.yaml (the set
                          of the mode) and, where the script
                          states the expected order and a series has three or more
                          refinements, the last rate with it; exit 1 on a mismatch, a missing
                          record or a run that raised (the ctests run --quick --check)
               --record   write the expected records of the mode from this run
The legacy positional `True` of the cython-era scripts means --quick.

The modules are built by CMake (reproducibles_python/CMakeLists.txt) into
<build>/reproducibles_python/; <build> is $HYPERHDG_BUILD_DIR or build/ next to the sources.
"""
import inspect
import math
import os
import random
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "python"))   # hyperhdg.py is a script next to the examples, not a package yet
import hyperhdg  # noqa: E402

# SciPy renamed the tolerance of the iterative solvers from tol to rtol in version 1.12.
_TOLERANCE = "rtol" if "rtol" in inspect.signature(spla.cg).parameters else "tol"

# the records written in this process, per script name: [(params, results)]
RECORDS = {}
# set by main() for --quick: the time-stepping scripts take far fewer steps then (see their test
# functions), so the expected files hold quick-mode records
QUICK = False
# the parameter that halves the mesh width from one record to the next
REFINEMENT_KEYS = ("refinement", "iteration")


def module(name):
  """The precompiled module reproducibles_python/modules/<name>.cxx, subscriptable by class name."""
  return hyperhdg.import_prebuilt(name, subdir="reproducibles_python")


# ---- linear algebra -------------------------------------------------------------------------

def operator(hdg, time=None):
  """The condensed system as a matrix-free scipy LinearOperator around trace_to_flux."""
  n = hdg.size_of_system()
  matvec = (lambda x: hdg.trace_to_flux(list(x))) if time is None \
      else (lambda x: hdg.trace_to_flux(list(x), time))
  return spla.LinearOperator((n, n), matvec=matvec)


def solve(A, rhs, tol=1e-13, solvers=("cg", "gmres", "bicgstab")):
  """Solve A x = rhs with the first solver; on failure fall back along the chain with the
  messages the cython-era scripts printed."""
  for k, name in enumerate(solvers):
    solution, info = getattr(spla, name)(A, rhs, **{_TOLERANCE: tol})
    if info == 0:
      return solution
    if k + 1 < len(solvers):
      print(f"{name.upper()} failed with a total number of {info} iterations. Trying {solvers[k + 1].upper()}!")
    else:
      print(f"{name.upper()} also failed with a total number of {info} iterations.")
  raise RuntimeError("Linear solvers did not converge!")


def time_stepping(hdg, time_steps, delta_time, tol=1e-13, solvers=("cg", "gmres", "bicgstab"),
                  set_final_data=True, on_start=None, on_step=None):
  """The time loop of the parabolic scripts: from make_initial, in every step the residual of the
  zero vector at the new time is the right-hand side, the condensed system is solved, and the
  solution becomes the data of the next step (set_final_data=False skips that after the last
  step, as the bilaplacian script did). on_start(solution) is called before the loop and
  on_step(step, solution) after every step (the movies plot there). Returns the solution at
  time_steps * delta_time."""
  solution = hdg.make_initial(hdg.zero_vector())
  if on_start:
    on_start(solution)
  A = operator(hdg)
  for step in range(time_steps):
    time = (step + 1) * delta_time
    rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector(), time))
    solution = solve(A, rhs, tol, solvers)
    if set_final_data or step + 1 < time_steps:
      hdg.set_data(list(solution), time)
    if on_step:
      on_step(step + 1, solution)
  return solution


# ---- eigenvalue problems ----------------------------------------------------------------------

class EigenvalueHelper:
  """The condensed operators of the eigenvalue scripts on the free dofs (the Dirichlet dofs,
  hdg.dirichlet_indices(), are removed from the vectors): stiffness = -trace_to_flux, mass either
  from trace_to_mass_flux (mass="flux", the MassApproxEigenvalue loop) or as the difference of
  trace_to_flux at the shift and at zero (mass="shift", the ShiftedInverseEigenvalue loop), and
  the shifted inverse (stiffness - sigma mass)^-1 by CG, for scipy's eigsh in shift-invert mode."""

  def __init__(self, hdg, sigma, mass):
    self.hdg, self.sigma, self.mass_kind = hdg, sigma, mass
    self.dirichlet = list(hdg.dirichlet_indices())
    n_long = hdg.size_of_system()
    self.index_vector = [-1] * (n_long - len(self.dirichlet))   # short index -> long index
    self.vector_index = [-1] * n_long                           # long index -> short index or -1
    n_indices = 0
    for i in range(len(self.index_vector)):
      while n_indices < len(self.dirichlet) and i + n_indices == self.dirichlet[n_indices]:
        n_indices += 1
      self.index_vector[i] = i + n_indices
      self.vector_index[i + n_indices] = i
    assert -1 not in self.index_vector

  def long_vector(self, vector):
    return [vector[x] if x > -1 else 0. for x in self.vector_index]

  def short_vector(self, vector):
    return [vector[x] for x in self.index_vector]

  def multiply_stiff(self, vector):
    return np.multiply(self.short_vector(self.hdg.trace_to_flux(self.long_vector(vector))), -1.)

  def multiply_mass(self, vector):
    if self.mass_kind == "flux":
      return self.short_vector(self.hdg.trace_to_mass_flux(self.long_vector(vector)))
    long = self.long_vector(vector)
    vec = np.subtract(self.hdg.trace_to_flux(long), self.hdg.trace_to_flux(long, self.sigma))
    return np.multiply(self.short_vector(vec), -1. / self.sigma)

  def shifted_mult(self, vector):
    if self.mass_kind == "flux":
      return np.subtract(self.multiply_stiff(vector), np.multiply(self.multiply_mass(vector), self.sigma))
    return np.multiply(self.short_vector(self.hdg.trace_to_flux(self.long_vector(vector), self.sigma)), -1.)

  def shifted_inverse(self, vector):
    if not hasattr(self, "shift_op"):
      n = len(self.index_vector)
      self.shift_op = spla.LinearOperator((n, n), matvec=self.shifted_mult)
    vec, info = spla.cg(self.shift_op, vector, **{_TOLERANCE: 1e-9})
    assert info == 0
    return vec

  def eigsh(self, k=1):
    """The k eigenvalues closest to the shift and their eigenvectors on all dofs."""
    n = len(self.index_vector)
    stiff = spla.LinearOperator((n, n), matvec=self.multiply_stiff)
    mass = spla.LinearOperator((n, n), matvec=self.multiply_mass)
    shifted_inv = spla.LinearOperator((n, n), matvec=self.shifted_inverse)
    vals, vecs = spla.eigsh(stiff, k=k, M=mass, sigma=self.sigma, which="LM", OPinv=shifted_inv)
    return vals, [self.long_vector([x[j].real for x in vecs]) for j in range(k)]


class NewtonHelper:
  """Residual and Jacobian of the nonlinear eigenvalue problem of the NonlinearEigenvalue loop:
  the unknown is (eigenvector, eigenvalue) in one vector, the last equation is the normalization
  |eigenvector|^2 = 1."""

  def __init__(self, hdg):
    self.hdg = hdg
    self.val = [0] * (hdg.size_of_system() + 1)

  def set_val(self, val):
    self.val = val

  def eval_residual(self, vector):
    vec = self.hdg.trace_to_flux(list(vector), vector[-1])
    temp, vector[-1] = vector[-1], 0.
    vec[-1] = np.linalg.norm(vector) ** 2 - 1.
    vector[-1] = temp
    return vec

  def eval_jacobian(self, vector):
    vec = self.hdg.jacobian_of_trace_to_flux(list(vector), vector[-1], list(self.val), self.val[-1])
    temp, self.val[-1] = self.val[-1], 0.
    vec[-1] = 2. * np.dot(vector, self.val)
    self.val[-1] = temp
    return vec


def newton_eigenvalue(hdg, initial, guess, n_newton_steps, solvers, tol_cap=None,
                      scaling_fac=1e-2, norm_exact=1e-9, alpha=0.1, beta=0.5):
  """Newton's method with backtracking line search on the nonlinear eigenvalue problem. initial
  is "default" (make_initial, normalized, eigenvalue guess with a random perturbation of up to
  1e-1 as the scripts always did) or a start vector with the eigenvalue as last entry. solvers:
  the Krylov solvers tried in turn for the Newton update, with tolerance scaling_fac * residual
  norm (capped at tol_cap if given). Returns the eigenvalue and the vector."""
  helper = NewtonHelper(hdg)
  n = hdg.size_of_system() + 1
  A = spla.LinearOperator((n, n), matvec=helper.eval_jacobian)
  if initial == "default":
    vector = np.asarray(hdg.make_initial([0] * n))
    vector /= np.linalg.norm(vector)
    vector[-1] = guess + 1e-3 * random.randint(-100, 100)
  else:
    vector = np.asarray(initial, dtype=float)
    temp, vector[-1] = vector[-1], 0.
    vector /= np.linalg.norm(vector)
    vector[-1] = temp
  residual = helper.eval_residual(vector)
  norm_res = np.linalg.norm(residual)
  for newton_step in range(n_newton_steps):
    helper.set_val(vector)
    norm_old, gamma = norm_res, 1
    tol = scaling_fac * norm_res if tol_cap is None else min(tol_cap, scaling_fac * norm_res)
    update = solve(A, residual, tol, solvers)
    candidate = np.subtract(vector, update)
    residual = helper.eval_residual(candidate)
    norm_res = np.linalg.norm(residual)
    while norm_res > (1 - alpha * gamma) * norm_old:
      if gamma < 1e-4:
        raise RuntimeError("Newton step is too small!")
      gamma = beta * gamma
      candidate = np.subtract(vector, np.multiply(update, gamma))
      residual = helper.eval_residual(candidate)
      norm_res = np.linalg.norm(residual)
    vector = candidate
    if norm_res < norm_exact:
      print("Newton solver converged after ", newton_step + 1, " iterations.")
      break
  if norm_res >= norm_exact:
    raise RuntimeError("Newton solver did not converge!")
  return vector[-1], vector


# ---- records ----------------------------------------------------------------------------------

def _value(x):
  """One YAML scalar: ints as they are, floats with all digits, everything else as str."""
  if isinstance(x, (bool, np.bool_)):
    return str(bool(x)).lower()
  if isinstance(x, (int, np.integer)):
    return str(int(x))
  if isinstance(x, (float, np.floating)):
    return repr(float(x))
  return str(x)


def _record(params, results):
  return "- {" + ", ".join(f"{k}: {_value(v)}" for k, v in {**params, **results}.items()) + "}"


def _parse(text):
  """The records of a YAML list of flow mappings as written by report(): [fields dict]."""
  records = []
  for line in text.splitlines():
    line = line.strip()
    if not line.startswith("- {"):
      continue
    fields = {}
    for key, value in re.findall(r"(\w+): ([^,}]+)", line[3:-1]):
      for convert in (int, float, str):
        try:
          fields[key] = convert(value)
          break
        except ValueError:
          pass
    records.append(fields)
  return records


def _normalize(params):
  """Parameters as they come back from the file, so records compare equal."""
  return _parse(_record(params, {}))[0]


def report(name, results, **params):
  """One run's record: results is the error (a number, key "error") or a dict of named results.
  Printed and appended to output/<name>.yaml; kept for the summary, --check and --record."""
  results = {"error": results} if not isinstance(results, dict) else results
  line = _record(params, results)
  print(line, flush=True)
  os.makedirs("output", exist_ok=True)
  with open(f"output/{name}.yaml", "a") as f:
    f.write(line + "\n")
  RECORDS.setdefault(name, []).append((_normalize(params), dict(results)))


def plot(hdg, solution, filename, time=None, options={}):
  """Plot a solution to output/<filename>.vtu with the options the reproducibles use, plus
  the script's own (e.g. the bilaplacian's boundaryScale and plotEdgeBoundaries)."""
  for option, value in {"fileName": filename, "printFileNumber": "false", "scale": "0.95",
                        **options}.items():
    hdg.plot_option(option, value)
  if time is None:
    hdg.plot_solution(list(solution))
  else:
    hdg.plot_solution(list(solution), time)


# ---- convergence rates ------------------------------------------------------------------------

def series(records):
  """Group records by every parameter but the refinement one: {(fixed params): [(level, results)]},
  the levels sorted; the refinement key, or None if there is none."""
  key = next((k for k in REFINEMENT_KEYS if any(k in p for p, _ in records)), None)
  if key is None:
    return {}, None
  groups = {}
  for params, results in records:
    fixed = tuple((k, v) for k, v in params.items() if k != key)
    groups.setdefault(fixed, []).append((params[key], results))
  return {k: sorted(v, key=lambda t: t[0]) for k, v in groups.items()}, key


def rates(levels):
  """The experimental orders log2(e_k / e_k+1) of each result along the sorted levels, per result
  key: {key: [rate, ...]} (None where an error is zero or missing)."""
  out = {}
  for (_, a), (_, b) in zip(levels, levels[1:]):
    for k in a:
      x, y = a.get(k), b.get(k)
      ok = isinstance(x, (int, float)) and isinstance(y, (int, float)) and x > 0 and y > 0
      out.setdefault(k, []).append(math.log2(x / y) if ok else None)
  return out


def summary(name, order=None):
  """Print the errors and rates of every series of <name>; return the complaints against the
  expected order (a number or a function of the fixed parameters) for series with three or
  more levels: the last rate must reach it up to ORDER_TOL."""
  groups, key = series(RECORDS.get(name, []))
  complaints = []
  if key is None:
    return complaints
  print(f"\nconvergence of {name} over {key}:")
  for fixed, levels in groups.items():
    params = dict(fixed)
    label = ", ".join(f"{k}={v}" for k, v in fixed)
    for rkey, r in rates(levels).items():
      errors = " ".join(f"{res.get(rkey):.3e}" if isinstance(res.get(rkey), (int, float)) else "-" for _, res in levels)
      rate_text = " ".join("-" if x is None else f"{x:.2f}" for x in r)
      expected = order(params) if callable(order) else order
      text = f"  {label} {rkey}: {errors} | rates {rate_text}"
      if expected is not None:
        text += f" | expected {expected:.2f}"
        if len(levels) >= 3 and (r[-1] is None or r[-1] < expected - ORDER_TOL):
          complaints.append(f"{name} {label} {rkey}: last rate {r[-1]} below expected {expected}")
          text += " FAILED"
      print(text)
  return complaints


ORDER_TOL = 0.25


# ---- expected records -------------------------------------------------------------------------

def expected_path(name):
  """expected/quick/<name>.yaml in quick mode, expected/full/<name>.yaml otherwise: the two
  sets differ where the time-stepping scripts take fewer steps in quick mode."""
  return HERE / "expected" / ("quick" if QUICK else "full") / f"{name}.yaml"


def check(name, rtol, atol):
  """Compare this run's records of <name> with the expected ones (expected_path): a record with
  the same parameters must exist and every result must agree to atol + rtol * |expected|.
  Returns the list of complaints (empty: all good)."""
  path = expected_path(name)
  if not path.exists():
    return [f"{path.relative_to(ROOT)} is missing (run with --record to write it)"]
  expected = _parse(path.read_text())
  complaints = []
  for params, results in RECORDS.get(name, []):
    match = [e for e in expected if all(e.get(k) == v for k, v in params.items())]
    if not match:
      complaints.append(f"{_record(params, results)}: no expected record")
      continue
    for key, value in results.items():
      ref = match[0].get(key)
      if not isinstance(ref, (int, float)) or abs(value - ref) > atol + rtol * abs(ref):
        complaints.append(f"{_record(params, results)}: expected {key}: {ref}")
  return complaints


def record(name):
  """Write the expected records of <name> (expected_path) from this run."""
  path = expected_path(name)
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text("".join(_record(p, r) + "\n" for p, r in RECORDS.get(name, [])))
  print(f"wrote {path.relative_to(ROOT)}: {len(RECORDS.get(name, []))} records")


# ---- the loops ------------------------------------------------------------------------------

def main(test, argv=None, rtol=1e-6, atol=1e-12, order=None, **grids):
  """Run test(**params) over the product of the grids (keyword order = loop nesting). --quick
  keeps only the first two values of the refinement grid (refinement or iteration, else the
  innermost), so every class runs once on coarse meshes. A RuntimeError of one run is printed as
  ERROR and the loops go on. Afterwards the convergence summary; order is the expected rate (a
  number or a function of the fixed parameters). --check compares the records with the expected
  files (tolerances rtol, atol) and the rates with order; --record writes the expected files;
  with --check the process exits 1 on any complaint or raised run."""
  global QUICK
  argv = sys.argv[1:] if argv is None else argv
  quick = QUICK = "--quick" in argv or "True" in argv
  checking, recording = "--check" in argv, "--record" in argv
  grids = {key: list(values) for key, values in grids.items()}
  if quick:
    key = next((k for k in REFINEMENT_KEYS if k in grids), next(reversed(grids)))
    grids[key] = grids[key][:2]
  failures = []

  def loops(keys, params):
    if not keys:
      start = datetime.now()
      print("Starting time is", start)
      try:
        test(**params)
      except RuntimeError as error:
        print("ERROR: ", error)
        failures.append(f"{params}: {error}")
      end = datetime.now()
      print("Program ended at", end, "after", end - start)
      return
    key = keys[0]
    for value in grids[key]:
      print(f"\n{key} = {value}\n")
      loops(keys[1:], {**params, key: value})

  loops(list(grids), {})
  complaints = list(failures)
  for name in RECORDS:
    complaints += summary(name, order)
  if recording:
    for name in RECORDS:
      record(name)
  if checking:
    for name in RECORDS:
      complaints += check(name, rtol, atol)
    for c in complaints:
      print("CHECK FAILED:", c)
    print(f"check: {len(complaints)} complaints in {sum(len(r) for r in RECORDS.values())} records")
    if complaints:
      sys.exit(1)
