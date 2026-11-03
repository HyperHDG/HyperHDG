#!/usr/bin/env python3
"""The smallest eigenvalue of the Laplacian (DiffusionEigs, NonlinearEigenvalue loop) on the cross
network domains/cross_2d.geo with every edge split into 2^iteration pieces, degrees 1 to 3,
iterations 2 to 5: Newton's method with line search on the eigenpair, at most 10 steps, GMRES then
BiCGStab for the updates, started from the loop's initial guess with dimension pi^2 (perturbed at
random) as the eigenvalue. The eigenvalues go to
output/diffusion_hypergraph_convergence_eigenvalue_newt_cross.yaml, the
eigenvectors are plotted. The classes DiffusionEigs_D1_P<degree> come from
modules/diffusion_file_nonlinear_eigenvalue.cxx.

    python3 reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_newt_cross.py [--quick]
"""
import numpy as np

import common

NAME = "diffusion_hypergraph_convergence_eigenvalue_newt_cross"
hdg_module = common.module("diffusion_file_nonlinear_eigenvalue")
DOMAIN = str(common.ROOT / "domains" / "cross_2d.geo")


def eigenvalue_newt(poly_degree, dimension, iteration, initial="default"):
  hdg = hdg_module[f"DiffusionEigs_D1_P{poly_degree}"](DOMAIN)
  hdg.set_refinement(2 ** iteration)
  value, vector = common.newton_eigenvalue(hdg, initial, dimension * (np.pi ** 2), n_newton_steps=10,
                                           solvers=("gmres", "bicgstab"))
  common.report(NAME, {"ev": value}, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, vector, f"diff_eig_newt-{dimension}-{iteration}", time=value)
  return value, vector, value


if __name__ == "__main__":
  common.main(eigenvalue_newt, rtol=1e-3, atol=1e-8, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(2, 6))
