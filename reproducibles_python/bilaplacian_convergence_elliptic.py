#!/usr/bin/env python3
"""Convergence of the LDG-H bilaplacian solver (Bilaplacian, Elliptic loop) on the unit cube in one
and two dimensions, degrees 1 to 3, six refinements, CG tolerance 1e-9: the errors go to
output/bilaplacian_convergence_elliptic.yaml, every run is plotted. The classes
Bilaplacian_D<dim>_P<degree> come from modules/bilaplacian_elliptic.cxx.

    python3 reproducibles_python/bilaplacian_convergence_elliptic.py [--quick]
"""
import numpy as np

import common

NAME = "bilaplacian_convergence_elliptic"
hdg_module = common.module("bilaplacian_elliptic")


def bilaplacian_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"Bilaplacian_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs, tol=1e-9)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"bil_conv_ellip-{dimension}-{iteration}", options={"boundaryScale": "0.9", "plotEdgeBoundaries": "true"})


if __name__ == "__main__":
  common.main(bilaplacian_test, rtol=1e-2, order=lambda p: p['poly_degree']  + 1, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
