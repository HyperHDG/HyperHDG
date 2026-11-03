#!/usr/bin/env python3
"""Convergence of the LDG-H bilaplacian solver (Bilaplacian, Elliptic loop) on hypergraphs of one-
and two-dimensional hyperedges in the unit cube, degrees 1 to 3, six refinements: the errors go to
output/bilaplacian_hypergraph_convergence_elliptic.yaml, every run is plotted. The classes
Bilaplacian_D<dim>_P<degree> come from modules/bilaplacian_hypergraph_elliptic.cxx.

    python3 reproducibles_python/bilaplacian_hypergraph_convergence_elliptic.py [--quick]
"""
import numpy as np

import common

NAME = "bilaplacian_hypergraph_convergence_elliptic"
hdg_module = common.module("bilaplacian_hypergraph_elliptic")


def bilaplacian_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"Bilaplacian_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"bil_conv_hyg-{dimension}-{iteration}", options={"boundaryScale": "0.9", "plotEdgeBoundaries": "true"})


if __name__ == "__main__":
  common.main(bilaplacian_test, rtol=1e-2, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
