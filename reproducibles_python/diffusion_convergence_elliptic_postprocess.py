#!/usr/bin/env python3
"""Convergence of the postprocessed LDG-H diffusion solver (DiffusionPostprocess, Elliptic loop) on
the unit cube in one and two dimensions, degrees 1 to 3, six refinements: the errors (the
postprocessed one first) go to output/diffusion_convergence_elliptic_postprocess.yaml, every run is
plotted. The classes DiffusionPostprocess_D<dim>_P<degree> come from
modules/diffusion_postprocess_elliptic.cxx.

    python3 reproducibles_python/diffusion_convergence_elliptic_postprocess.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_elliptic_postprocess"
hdg_module = common.module("diffusion_postprocess_elliptic")


def diffusion_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"DiffusionPostprocess_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_conv_ellip-{dimension}-{iteration}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, order=lambda p: p['poly_degree']  + 2, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
