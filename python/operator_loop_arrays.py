#!/usr/bin/env python3
"""Check that the in-place methods of the local-operator loop reject output arrays they could
not write into.

    python python/operator_loop_arrays.py

apply, residual and assemble write into caller-owned NumPy arrays. An output array of the wrong
dtype or a non-contiguous view would be converted to a temporary by nanobind, and the result
silently lost; the bindings therefore declare these arguments noconvert, so such calls raise
TypeError. Inputs may still be converted.
"""

import sys

import numpy as np

import operator_loop


def raises(call):
  try:
    call()
  except TypeError:
    return True
  return False


def main():
  op = operator_loop.problem(operator_loop.Options(dim=1, degree=1, n=8))
  n, m = op.size_of_system(), op.n_triplets(True)
  x, y = np.ones(n), np.empty(n)
  rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)

  op.apply(x, y)
  op.residual(x, y)
  op.assemble(rows, cols, vals, 0., True)
  op.apply(np.ones(n, np.float32), y)  # inputs are converted

  checks = {
    "apply wrong dtype": lambda: op.apply(x, np.empty(n, np.float32)),
    "apply strided": lambda: op.apply(x, np.empty(2 * n)[::2]),
    "apply list": lambda: op.apply(x, [0.] * n),
    "residual wrong dtype": lambda: op.residual(x, np.empty(n, np.float32)),
    "assemble wrong index dtype": lambda: op.assemble(np.empty(m, np.int64), cols, vals, 0., True),
    "assemble wrong value dtype": lambda: op.assemble(rows, cols, np.empty(m, np.float32), 0., True),
  }
  failed = [name for name, call in checks.items() if not raises(call)]
  for name in failed:
    print("no TypeError:", name)
  print("ok" if not failed else "failed")
  return 1 if failed else 0


if __name__ == "__main__":
  sys.exit(main())
