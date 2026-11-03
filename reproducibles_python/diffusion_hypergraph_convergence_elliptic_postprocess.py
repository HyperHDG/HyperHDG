#!/usr/bin/env python3
"""Convergence of the postprocessed LDG-H diffusion solver (DiffusionPostprocess, Elliptic loop) on
hypergraphs of one- to three-dimensional hyperedges in the unit cube, degrees 1 to 3, six
refinements: the errors go to output/diffusion_hypergraph_convergence_elliptic_postprocess.yaml,
every run is plotted. The classes DiffusionPostprocess_D<dim>_P<degree> come from
modules/diffusion_postprocess_hypergraph_elliptic.cxx.

    python3 reproducibles_python/diffusion_hypergraph_convergence_elliptic_postprocess.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_hypergraph_convergence_elliptic_postprocess"
hdg_module = common.module("diffusion_postprocess_hypergraph_elliptic")


def diffusion_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"DiffusionPostprocess_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_conv_hyg-{dimension}-{iteration}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, poly_degree=range(1, 4), dimension=range(1, 4), iteration=range(6))
