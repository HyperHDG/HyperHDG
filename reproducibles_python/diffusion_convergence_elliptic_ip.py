#!/usr/bin/env python3
"""Convergence of the interior penalty diffusion solver (DiffusionIP, Elliptic loop) on the unit
cube in one and two dimensions, degrees 1 to 3, six refinements, penalty 5 (d+1) p^3 2^it: the
errors go to output/diffusion_convergence_elliptic_ip.yaml, every run is plotted. The classes
DiffusionIP_D<dim>_P<degree> come from modules/diffusion_ip_elliptic.cxx.

    python3 reproducibles_python/diffusion_convergence_elliptic_ip.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_elliptic_ip"
hdg_module = common.module("diffusion_ip_elliptic")


def diffusion_test(poly_degree, dimension, iteration):
  penalty = 5 * (dimension + 1) * (poly_degree ** 3) * (2 ** iteration)
  hdg = hdg_module[f"DiffusionIP_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension, penalty)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_conv_ellip-{dimension}-{iteration}")


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-3, order=lambda p: p['poly_degree']  + 1, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
