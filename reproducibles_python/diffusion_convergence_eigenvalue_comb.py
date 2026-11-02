#!/usr/bin/env python3
"""Newton's method on the eigenpair started from the shifted-inverse (SI) and from the
mass-approximation (MA) eigenpair of the same discretization: the errors of the start and of the
Newton result go to output/diffusion_convergence_eigenvalue_comb.yaml.

    python3 reproducibles_python/diffusion_convergence_eigenvalue_comb.py [--quick]
"""
from diffusion_convergence_eigenvalue_approx import eigenvalue_approx_MA
from diffusion_convergence_eigenvalue_shifted_inverse import eigenvalue_approx_SI
from diffusion_convergence_eigenvalue_newt import eigenvalue_newt

import common

NAME = "diffusion_convergence_eigenvalue_comb"


def combination(poly_degree, dimension, iteration, approx_type):
  start, name = (eigenvalue_approx_SI, "SI") if approx_type == 0 else (eigenvalue_approx_MA, "MA")
  value, vector, error_start = start(poly_degree, dimension, iteration)
  vector.append(value)
  try:
    _, _, error_newton = eigenvalue_newt(poly_degree, dimension, iteration, vector)
  except RuntimeError as error:
    print("ERROR: ", error)
    error_newton = -1.
  common.report(NAME, {"error_start": error_start, "error_newton": error_newton},
                poly_degree=poly_degree, dimension=dimension, iteration=iteration, start=name)


if __name__ == "__main__":
  common.main(combination, rtol=1e-3, atol=1e-8, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6),
              approx_type=range(2))
