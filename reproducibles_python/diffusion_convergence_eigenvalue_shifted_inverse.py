#!/usr/bin/env python3
"""The smallest eigenvalue of the Laplacian (DiffusionEigs, ShiftedInverseEigenvalue loop) on the unit cube in one and two dimensions,
degrees 1 to 3, refinements 2 to 5: scipy's eigsh in shift-invert mode around half the exact eigenvalue, the mass
matrix recovered from the loop's shifted operator. The eigenvalue errors go to
output/diffusion_convergence_eigenvalue_shifted_inverse.yaml, the eigenvectors are plotted. The classes DiffusionEigs_D<dim>_P<degree> come from
modules/diffusion_shifted_inverse_eigenvalue.cxx.

    python3 reproducibles_python/diffusion_convergence_eigenvalue_shifted_inverse.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_eigenvalue_shifted_inverse"
hdg_module = common.module("diffusion_shifted_inverse_eigenvalue")


def eigenvalue_approx_SI(poly_degree, dimension, iteration):
  exact_eigenval = dimension * (np.pi ** 2)
  sigma = exact_eigenval / 2.
  hdg = hdg_module[f"DiffusionEigs_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  vals, vecs = common.EigenvalueHelper(hdg, sigma, mass="shift").eigsh()
  error = np.absolute(vals[0] - exact_eigenval)
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vecs[0], f"diff_eig_shifi-{dimension}-{iteration}", time=vals[0].real)
  return vals[0].real, vecs[0], error


if __name__ == "__main__":
  common.main(eigenvalue_approx_SI, rtol=1e-3, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
