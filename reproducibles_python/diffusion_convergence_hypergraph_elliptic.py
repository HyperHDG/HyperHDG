#!/usr/bin/env python3
"""Convergence of the LDG-H diffusion solver (Diffusion, Elliptic loop) on hypergraphs of one- to
three-dimensional hyperedges in the unit cube under refinement of the hyperedges (2^refinement
subintervals per hyperedge on a 4 x 4 x 4 cube), degrees 1 to 3: the errors go to
output/diffusion_convergence_hypergraph_elliptic.yaml, every run is plotted. The classes
Diffusion_D<dim>_P<degree> come from modules/diffusion_hypergraph_elliptic.cxx.

    python3 reproducibles_python/diffusion_convergence_hypergraph_elliptic.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_hypergraph_elliptic"
hdg_module = common.module("diffusion_hypergraph_elliptic")


def diffusion_test(poly_degree, dimension, iteration, refinement):
  hdg = hdg_module[f"Diffusion_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  hdg.set_refinement(2 ** refinement)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration, refinement=refinement)
  common.plot(hdg, solution, f"diff_hyg_conv-{dimension}-{refinement}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, poly_degree=range(1, 4), dimension=range(1, 4), iteration=[2], refinement=range(4))
