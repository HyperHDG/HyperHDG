#!/usr/bin/env python3
"""Convergence of the Timoshenko beam solver (TimoshenkoBeam with its default parameters, Elliptic
loop) on the network domains/cross.geo with every edge split into 2^iteration pieces, degrees 1 to
6, iterations 0 to 5, the assembled system solved by CG to 1e-12: the errors go to
output/timoshenko_convergence_cross.yaml, every run is plotted. The classes TimoshenkoBeam_D1_P<degree>
come from modules/timoshenko_file_elliptic.cxx.

    python3 reproducibles_python/timoshenko_convergence_cross.py [--quick]
"""
import numpy as np
import scipy.sparse as sp

import common

NAME = "timoshenko_convergence_cross"
hdg_module = common.module("timoshenko_file_elliptic")
DOMAIN = str(common.ROOT / "domains" / "cross.geo")


def diffusion_test(poly_degree, iteration):
  hdg = hdg_module[f"TimoshenkoBeam_D1_P{poly_degree}"](DOMAIN)
  hdg.set_refinement(2 ** iteration)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  n = hdg.size_of_system()
  rows, cols, vals = hdg.trace_to_flux_mat()
  A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
  solution = common.solve(A, rhs, tol=1e-12, solvers=("cg",))
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, iteration=iteration)
  common.plot(hdg, solution, f"timo_cross-{poly_degree}-{iteration}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, poly_degree=range(1, 7), iteration=range(6))
