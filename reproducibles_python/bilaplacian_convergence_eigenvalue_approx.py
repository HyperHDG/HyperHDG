#!/usr/bin/env python3
"""The smallest eigenvalue of the bilaplacian (Bilaplacian, MassApproxEigenvalue loop) on the unit cube in one and two dimensions,
degrees 1 to 3, refinements 2 to 5: scipy's eigsh in shift-invert mode around half the exact
eigenvalue (three quarters for the bilaplacian), the mass from the loop's trace_to_mass_flux. The eigenvalue errors go to
output/bilaplacian_convergence_eigenvalue_approx.yaml, the eigenvectors are plotted. The classes Bilaplacian_D<dim>_P<degree> come from
modules/bilaplacian_mass_approx_eigenvalue.cxx.

    python3 reproducibles_python/bilaplacian_convergence_eigenvalue_approx.py [--quick]
"""
import numpy as np

import common

NAME = "bilaplacian_convergence_eigenvalue_approx"
hdg_module = common.module("bilaplacian_mass_approx_eigenvalue")


def eigenvalue_approx_MA(poly_degree, dimension, iteration):
  exact_eigenval = (dimension * (np.pi ** 2)) ** 2
  sigma = 3. * exact_eigenval / 4.
  hdg = hdg_module[f"Bilaplacian_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  vals, vecs = common.EigenvalueHelper(hdg, sigma, mass="flux").eigsh()
  error = np.absolute(vals[0] - exact_eigenval)
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vecs[0], f"bil_e-{dimension}-{iteration}")
  return vals[0].real, vecs[0], error


if __name__ == "__main__":
  common.main(eigenvalue_approx_MA, rtol=1e-3, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6) if not common.QUICK else [2])
