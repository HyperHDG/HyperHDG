#!/usr/bin/env python3
"""The smallest eigenvalue of the Laplacian (Diffusion, MassApproxEigenvalue loop) on hypergraphs of one- and two-dimensional hyperedges in the unit cube,
degrees 1 to 3, refinements 2 to 5: scipy's eigsh in shift-invert mode around half the exact
eigenvalue, the mass from the loop's trace_to_mass_flux. The eigenvalue errors go to
output/diffusion_hypergraph_convergence_eigenvalue_approx.yaml, the eigenvectors are plotted. The classes Diffusion_D<dim>_P<degree> come from
modules/diffusion_hypergraph_mass_approx_eigenvalue.cxx.

    python3 reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_approx.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_hypergraph_convergence_eigenvalue_approx"
hdg_module = common.module("diffusion_hypergraph_mass_approx_eigenvalue")


def eigenvalue_approx_MA(poly_degree, dimension, iteration):
  exact_eigenval = dimension * (np.pi ** 2)
  sigma = exact_eigenval / 2.
  hdg = hdg_module[f"Diffusion_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  vals, vecs = common.EigenvalueHelper(hdg, sigma, mass="flux").eigsh()
  error = np.absolute(vals[0] - exact_eigenval)
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vecs[0], f"diff_e-{dimension}-{iteration}")
  return vals[0].real, vecs[0], error


if __name__ == "__main__":
  common.main(eigenvalue_approx_MA, rtol=1e-3, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
