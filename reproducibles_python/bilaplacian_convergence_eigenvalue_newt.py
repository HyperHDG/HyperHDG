#!/usr/bin/env python3
"""The smallest eigenvalue of the bilaplacian (BilaplacianEigs, NonlinearEigenvalue loop) on the unit cube in one and two dimensions, degrees
1 to 3, refinements 2 to 5: Newton's method with line search on the eigenpair, at most
25 steps, BiCGStab then GMRES, tolerance capped at 1e-9, for the updates, started from the loop's initial guess with the
eigenvalue perturbed at random. The eigenvalue errors go to output/bilaplacian_convergence_eigenvalue_newt.yaml, the
eigenvectors are plotted. The classes BilaplacianEigs_D<dim>_P<degree> come from modules/bilaplacian_nonlinear_eigenvalue.cxx.

    python3 reproducibles_python/bilaplacian_convergence_eigenvalue_newt.py [--quick]
"""
import numpy as np

import common

NAME = "bilaplacian_convergence_eigenvalue_newt"
hdg_module = common.module("bilaplacian_nonlinear_eigenvalue")


def eigenvalue_newt(poly_degree, dimension, iteration, initial="default"):
  exact_eigenval = (dimension * (np.pi ** 2)) ** 2
  hdg = hdg_module[f"BilaplacianEigs_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension)
  value, vector = common.newton_eigenvalue(hdg, initial, exact_eigenval, n_newton_steps=25,
                                           solvers=("bicgstab", "gmres"), tol_cap=1e-9)
  error = np.absolute(value - exact_eigenval)
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vector, f"bil_eig_newt-{dimension}-{iteration}", time=value)
  return value, vector, error


if __name__ == "__main__":
  common.main(eigenvalue_newt, rtol=1e-3, atol=1e-8, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
