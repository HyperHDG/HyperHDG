#!/usr/bin/env python3
"""The smallest eigenvalue of the Laplacian (DiffusionEigs, NonlinearEigenvalue loop) on hypergraphs of one- and two-dimensional hyperedges in the unit cube, degrees
1 to 3, refinements 2 to 5: Newton's method with line search on the eigenpair, at most
10 steps, GMRES then BiCGStab for the updates, started from the loop's initial guess with the
eigenvalue perturbed at random. The eigenvalue errors go to output/diffusion_hypergraph_convergence_eigenvalue_newt.yaml, the
eigenvectors are plotted. The classes DiffusionEigs_D<dim>_P<degree> come from modules/diffusion_hypergraph_nonlinear_eigenvalue.cxx.

    python3 reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_newt.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_hypergraph_convergence_eigenvalue_newt"
hdg_module = common.module("diffusion_hypergraph_nonlinear_eigenvalue")


def eigenvalue_newt(poly_degree, dimension, iteration, initial="default"):
  exact_eigenval = dimension * (np.pi ** 2)
  hdg = hdg_module[f"DiffusionEigs_D{dimension}_P{poly_degree}"]([2 ** iteration] * 3)
  value, vector = common.newton_eigenvalue(hdg, initial, exact_eigenval, n_newton_steps=10,
                                           solvers=("gmres", "bicgstab"))
  error = np.absolute(value - exact_eigenval)
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vector, f"diff_eig_newt-{dimension}-{iteration}", time=value)
  return value, vector, error


if __name__ == "__main__":
  common.main(eigenvalue_newt, order=lambda p: 2 * p['poly_degree'] + 1, rtol=1e-3, atol=1e-8, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
