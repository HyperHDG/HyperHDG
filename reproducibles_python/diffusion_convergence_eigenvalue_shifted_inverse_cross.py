#!/usr/bin/env python3
"""The two smallest eigenvalues of the Laplacian (DiffusionEigs, ShiftedInverseEigenvalue loop) on
the cross network domains/cross_2d.geo with every edge split into 2^iteration pieces, degrees 1 to
3, iterations 2 to 5: scipy's eigsh in shift-invert mode around dimension pi^2 / 2 (the "dimension"
only sets the shift here). The eigenvalues go to
output/diffusion_convergence_eigenvalue_shifted_inverse_cross.yaml, the second eigenvector is plotted. The classes DiffusionEigs_D1_P<degree>
come from modules/diffusion_file_shifted_inverse_eigenvalue.cxx.

    python3 reproducibles_python/diffusion_convergence_eigenvalue_shifted_inverse_cross.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_convergence_eigenvalue_shifted_inverse_cross"
hdg_module = common.module("diffusion_file_shifted_inverse_eigenvalue")
DOMAIN = str(common.ROOT / "domains" / "cross_2d.geo")


def eigenvalue_approx_SI(poly_degree, dimension, iteration):
  sigma = dimension * (np.pi ** 2) / 2.
  hdg = hdg_module[f"DiffusionEigs_D1_P{poly_degree}"](DOMAIN)
  hdg.set_refinement(2 ** iteration)
  vals, vecs = common.EigenvalueHelper(hdg, sigma, mass="shift").eigsh(k=2)
  error = vals[0]
  common.report(NAME, {"ev1": vals[0], "ev2": vals[1]}, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vecs[1], f"diff_eig_shifi-{dimension}-{iteration}", time=vals[1].real)
  return vals[0].real, vecs[1], error


if __name__ == "__main__":
  common.main(eigenvalue_approx_SI, rtol=1e-3, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
