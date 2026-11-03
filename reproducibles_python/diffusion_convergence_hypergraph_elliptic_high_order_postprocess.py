#!/usr/bin/env python3
"""Convergence of the postprocessed LDG-H diffusion solver (DiffusionPostprocess, Elliptic loop) on
the hypergraph of two-dimensional hyperedges in the unit cube (2 x 2 x 2) under refinement of the
hyperedges, degrees 0 to 3, the high-order test solution: the errors go to
output/diffusion_convergence_hypergraph_elliptic_high_order_postprocess.yaml,
every run is plotted. The classes DiffusionPostprocess_D2_P<degree> come from
modules/diffusion_postprocess_hypergraph_elliptic_ho.cxx.

    python3 reproducibles_python/diffusion_convergence_hypergraph_elliptic_high_order_postprocess.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_hypergraph_elliptic_high_order_postprocess"
hdg_module = common.module("diffusion_postprocess_hypergraph_elliptic_ho")


def diffusion_test(poly_degree, dimension, iteration, refinement):
  hdg = hdg_module[f"DiffusionPostprocess_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  hdg.set_refinement(2 ** refinement)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration, refinement=refinement)
  common.plot(hdg, solution, f"diff_hg_conv_ho_pp-{dimension}-{refinement}")


if __name__ == "__main__":
  common.main(diffusion_test, order=lambda p: p['poly_degree'] + 2, rtol=1e-5, poly_degree=range(4), dimension=[2], iteration=[1], refinement=range(7))
