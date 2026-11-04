#!/usr/bin/env python3
"""The NumPy bindings of GlobalLoop::Generic (include/HyperHDG/bind_python.hxx).

On the module tests/bindings_module.cxx, which CMake builds (PYTHONPATH points at it):
apply, residual and assemble take NumPy arrays, and the assembled matrix is the matrix-free
operator. Their output arrays are noconvert: a wrong dtype, a strided view or a list raises
TypeError instead of a silently lost result. Inputs are converted.
"""

import sys

import numpy as np

import hyperhdg_bindings_test

failures = 0


def check(ok, what):
  global failures
  print(f"{'ok  ' if ok else 'FAIL'}  {what}")
  failures += not ok


hdg = hyperhdg_bindings_test.Diffusion([8])
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
