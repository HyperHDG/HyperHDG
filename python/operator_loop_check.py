#!/usr/bin/env python3
"""Check the local-operator loop against GlobalLoop::Elliptic to roundoff.

    python python/operator_loop_check.py [--dim 2 --degree 2 --n 4]

Compares, on the same mesh and solver parameters: the matrix-free product, the residual, the
assembled matrix, the prepared store against streaming, the CG solution and the L2 error. Prints
the deviations as YAML; exit status 1 if any exceeds its tolerance.
"""

import argparse
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import yaml

import operator_loop

CASES = [(1, 1, 8), (1, 3, 8), (2, 1, 4), (2, 2, 4), (2, 3, 3), (3, 1, 2), (3, 2, 2)]


def csr(triplets, n):
  rows, cols, vals = triplets
  return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def cg(matrix, rhs):
  x, info = spla.cg(matrix, rhs, rtol=1e-12, maxiter=100 * rhs.size)
  assert info == 0
  return x


def check(dim, degree, n_elements):
  opt = operator_loop.Options(dim=dim, degree=degree, n=n_elements)
  ell = operator_loop.problem(operator_loop.Options(**{**opt.__dict__, "loop": "elliptic"}))
  op = operator_loop.problem(opt)
  n = ell.size_of_system()
  x = np.random.default_rng(1).standard_normal(n)
  free = np.ones(n, dtype=bool)
  free[op.constraints()] = False
  zero = np.zeros(n)

  def op_apply(v):
    op.apply(v, y := np.empty(n))
    return y
  def op_residual(v):
    op.residual(v, r := np.empty(n))
    return r
  def op_assemble(identity_rows=False):
    m = op.n_triplets(identity_rows)
    t = (np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m))
    op.assemble(*t, 0., identity_rows)
    return csr(t, n)

  y_ell, y_op = ell.trace_to_flux(x), op_apply(x)
  scale = np.abs(y_ell).max()
  r_ell, r_op = ell.residual_flux(x), op_residual(x)
  a_ell, a_op = csr(ell.trace_to_flux_mat(), n), op_assemble()
  op.prepare()
  y_prep, a_prep = op_apply(x), op_assemble()
  sol_ell = cg(a_ell, -ell.residual_flux(zero))
  sol_op = cg(op_assemble(True), -op_residual(zero))
  err_ell, err_op = ell.errors(sol_ell)[0], op.errors(sol_op)[0]

  deviations = {  # name: (value, tolerance)
    "apply_vs_trace_to_flux": (np.abs(y_ell - y_op)[free].max() / scale, 1e-13),
    "apply_constrained_rows_identity": (np.abs(y_op[~free] - x[~free]).max(), 0.),
    "residual_vs_residual_flux": (np.abs(r_ell - r_op)[free].max() / np.abs(r_ell).max(), 1e-13),
    "assemble_vs_trace_to_flux_mat": (abs(a_ell - a_op).max() / abs(a_ell).max(), 1e-13),
    "prepared_vs_streaming_apply": (np.abs(y_prep - y_op).max() / scale, 1e-13),
    "prepared_vs_streaming_assemble": (abs(a_prep - a_op).max() / abs(a_ell).max(), 1e-13),
    "cg_solution": (np.abs(sol_ell - sol_op).max() / np.abs(sol_ell).max(), 1e-9),
    "l2_error": (abs(err_ell - err_op) / err_ell, 1e-9),
  }
  return {"dim": dim, "degree": degree, "n": n_elements, "dofs": n,
          "operators": op.n_operators(), "error": float(err_op),
          "ok": all(v <= tol for v, tol in deviations.values()),
          "deviations": {k: float(v) for k, (v, _) in deviations.items()}}


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--dim", type=int)
  parser.add_argument("--degree", type=int)
  parser.add_argument("--n", type=int)
  opt = parser.parse_args()
  cases = [(opt.dim, opt.degree, opt.n)] if opt.dim else CASES
  results = [check(*case) for case in cases]
  print(yaml.safe_dump(results, sort_keys=False), end="")
  return 0 if all(r["ok"] for r in results) else 1


if __name__ == "__main__":
  sys.exit(main())
