#!/usr/bin/env python3
"""Convergence of the LDG-H diffusion solver (Diffusion, Elliptic loop) on the unit cube in one
and two dimensions, degrees 1 to 3, six refinements: the errors go to
output/diffusion_convergence_elliptic.yaml, every run is plotted. The classes
Diffusion_D<dim>_P<degree> come from modules/diffusion_elliptic.cxx.

    python3 reproducibles_python/diffusion_convergence_elliptic.py [--quick]
"""
import numpy as np
import scipy.sparse as sp

import common

NAME = "diffusion_convergence_elliptic"
hdg_module = common.module("diffusion_elliptic")


def diffusion_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"Diffusion_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  n = hdg.size_of_system()
  rows, cols, vals = hdg.trace_to_flux_mat()
  A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
  solution = common.solve(A, rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_conv_ellip-{dimension}-{iteration}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, order=lambda p: p['poly_degree']  + 1, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
