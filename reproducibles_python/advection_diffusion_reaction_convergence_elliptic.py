#!/usr/bin/env python3
"""Convergence of the LDG diffusion-advection-reaction solver (DiffusionAdvectionReaction, Elliptic
loop) on the unit cube in one and two dimensions, degrees 1 to 3, six refinements: the errors go
to output/advection_diffusion_reaction_convergence_elliptic.yaml, every run is plotted. The classes DiffusionAdvectionReaction_D<dim>_P<degree> come
from modules/diffusion_advection_reaction_elliptic.cxx.

    python3 reproducibles_python/advection_diffusion_reaction_convergence_elliptic.py [--quick]
"""
import numpy as np

import common

NAME = "advection_diffusion_reaction_convergence_elliptic"
hdg_module = common.module("diffusion_advection_reaction_elliptic")


def adv_dif_reac_test(poly_degree, dimension, iteration):
  hdg = hdg_module[f"DiffusionAdvectionReaction_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution = common.solve(common.operator(hdg), rhs)
  error = hdg.errors(list(solution))[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_adv_reac_conv_ellip-{dimension}-{iteration}")


if __name__ == "__main__":
  common.main(adv_dif_reac_test, rtol=1e-5, order=lambda p: p['poly_degree']  + 1, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
