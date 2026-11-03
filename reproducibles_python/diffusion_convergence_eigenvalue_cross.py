#!/usr/bin/env python3
"""Newton's method on the eigenpair of the cross network started from the shifted-inverse
eigenpair of the same discretization: the eigenvalue of the start and of the Newton result go to
output/diffusion_convergence_eigenvalue_cross.yaml.

    python3 reproducibles_python/diffusion_convergence_eigenvalue_cross.py [--quick]
"""
from diffusion_convergence_eigenvalue_shifted_inverse_cross import eigenvalue_approx_SI
from diffusion_hypergraph_convergence_eigenvalue_newt_cross import eigenvalue_newt

import common

NAME = "diffusion_convergence_eigenvalue_cross"


def combination(poly_degree, dimension, iteration, approx_type):
  value, vector, error_start = eigenvalue_approx_SI(poly_degree, dimension, iteration)
  vector.append(value)
  try:
    _, _, error_newton = eigenvalue_newt(poly_degree, dimension, iteration, vector)
  except RuntimeError as error:
    print("ERROR: ", error)
    error_newton = -1.
  common.report(NAME, {"error_start": error_start, "error_newton": error_newton},
                poly_degree=poly_degree, dimension=dimension, iteration=iteration, start="SI")


if __name__ == "__main__":
  common.main(combination, rtol=1e-3, atol=1e-8, poly_degree=range(1, 4), dimension=[1], iteration=range(2, 6), approx_type=[0])
