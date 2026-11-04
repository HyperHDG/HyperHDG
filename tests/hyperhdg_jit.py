#!/usr/bin/env python3
"""python/hyperhdg.py and the bindings of GlobalLoop::Generic (include/HyperHDG/bind_python.hxx).

- A module source is compiled for the running interpreter, and the extension module is found and
  imported; its classes are reached by attribute and by subscript.
- Variants of a module (hash_name=True) get build directories of their own, so switching back to
  one does not rebuild it.
- apply, residual and assemble take NumPy arrays: the assembled matrix is the matrix-free
  operator. Their output arrays are noconvert: a wrong dtype, a strided view or a list raises
  TypeError instead of a silently lost result. Inputs are converted.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import hyperhdg  # noqa: E402

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>

using Diffusion = GlobalLoop::Generic<
  Topology::Cubic<1, 1>, Geometry::UnitCube<1, 1, double>, NodeDescriptor::Cubic<1, 1>,
  LocalSolver::Diffusion<1, {degree}, {quadrature}, tests::DiffusionSin, double> >;

NB_MODULE(hyperhdg_jit_test, m)
{{
  HyperHDG::bind_python<Diffusion, std::vector<unsigned int> >(m, "Diffusion");
}}
"""

failures = 0


def check(ok, what):
  global failures
  print(f"{'ok  ' if ok else 'FAIL'}  {what}")
  failures += not ok


def code(degree):
  return MODULE.format(degree=degree, quadrature=2 * degree)


# ---- compiling, finding and importing
module = hyperhdg.load(code(1))
check(module.Diffusion is module["Diffusion"], "the class by attribute and by subscript")

# ---- variants: separate build directories, no rebuild when switching back
paths = [Path(hyperhdg.compile(code(degree), hash_name=True)) for degree in (1, 2)]
mtime = paths[0].stat().st_mtime_ns
again = Path(hyperhdg.compile(code(1), hash_name=True))
check(paths[0] != paths[1] and paths[0].parent != paths[1].parent,
      "two variants are built in separate directories")
check(again == paths[0] and again.stat().st_mtime_ns == mtime,
      "switching back to a variant does not rebuild it")

# ---- the bindings of Generic: NumPy in place
hdg = module.Diffusion([8])
hdg.constrain_dirichlet_nodes()
hdg.prepare()
n = hdg.size_of_system()
x = np.random.default_rng(1).standard_normal(n)
y = np.empty(n)
hdg.apply(x, y)
m = hdg.n_triplets(True)
rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
hdg.assemble(rows, cols, vals, 0., True)
matrix = np.zeros((n, n))
np.add.at(matrix, (rows, cols), vals)
check(np.abs(matrix @ x - y).max() <= 1e-12 * np.abs(y).max(),
      "assemble into NumPy arrays = apply into a NumPy array")
r = np.empty(n)
hdg.residual(x, r)
check(np.isfinite(r).all(), "residual into a NumPy array")
y32 = np.empty(n)
hdg.apply(x.astype(np.float32), y32)
check(np.abs(y32 - matrix @ x.astype(np.float32)).max() <= 1e-6 * np.abs(y32).max(),
      "an input of another dtype is converted")


def raises_type_error(call):
  try:
    call()
  except TypeError:
    return True
  return False


for what, call in {
    "apply: output of the wrong dtype": lambda: hdg.apply(x, np.empty(n, np.float32)),
    "apply: strided output": lambda: hdg.apply(x, np.empty(2 * n)[::2]),
    "apply: a list as output": lambda: hdg.apply(x, [0.] * n),
    "residual: output of the wrong dtype": lambda: hdg.residual(x, np.empty(n, np.float32)),
    "assemble: indices of the wrong dtype": lambda: hdg.assemble(np.empty(m, np.int64), cols,
                                                                 vals, 0., True),
    "assemble: values of the wrong dtype": lambda: hdg.assemble(rows, cols,
                                                                np.empty(m, np.float32), 0., True),
}.items():
  check(raises_type_error(call), what + " raises TypeError")

sys.exit(1 if failures else 0)
