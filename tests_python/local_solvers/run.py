# Compile and run every local solver (with its global loop) once on a very coarse mesh.
#
# The cases reuse the test functions of the scripts in reproducibles_python with small parameters.
# This ensures that all local solvers compile with the used compiler (not only the ones used by the
# other tests) and that they can be run. Usage: python3 run.py <case name>

import importlib, os, sys

main_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + "/../..")
sys.path.append(main_dir + "/import")
sys.path.append(main_dir + "/reproducibles_python")
os.chdir(main_dir)

## Case name -> (reproducible module, test function, arguments).
CASES = {
  "Diffusion_Elliptic":
    ("diffusion_convergence_elliptic", "diffusion_test", (1, 2, 1)),
  "Diffusion_Elliptic_hypergraph":
    ("diffusion_hypergraph_convergence_elliptic", "diffusion_test", (1, 2, 1)),
  "DiffusionIP_Elliptic":
    ("diffusion_convergence_elliptic_ip", "diffusion_test", (1, 2, 1)),
  "DiffusionPostprocess_Elliptic":
    ("diffusion_convergence_elliptic_postprocess", "diffusion_test", (1, 2, 1)),
  "DiffusionAdvectionReaction_Elliptic":
    ("advection_diffusion_reaction_convergence_elliptic", "adv_dif_reac_test", (1, 2, 1)),
  "Bilaplacian_Elliptic":
    ("bilaplacian_convergence_elliptic", "bilaplacian_test", (1, 2, 1)),
  "TimoshenkoBeam_Elliptic":
    ("timoshenko_convergence_cross", "diffusion_test", (1, 0)),
  "DiffusionParab_Parabolic":
    ("diffusion_convergence_parabolic", "diffusion_test", (1, 1, 1)),
  "BilaplacianParab_Parabolic":
    ("bilaplacian_convergence_parabolic", "bilaplacian_test", (1, 1, 1)),
  "AdvectionParab_Parabolic":
    ("advection_convergence_parabolic", "diffusion_test", (1., 1, 1, 1)),
  "DiffusionEigs_NonlinearEigenvalue":
    ("diffusion_convergence_eigenvalue_newt", "eigenvalue_newt", (1, 1, 2)),
  "BilaplacianEigs_NonlinearEigenvalue":
    ("bilaplacian_convergence_eigenvalue_newt", "eigenvalue_newt", (1, 1, 2)),
  "DiffusionEigs_ShiftedInverseEigenvalue":
    ("diffusion_convergence_eigenvalue_shifted_inverse", "eigenvalue_approx_SI", (1, 1, 2)),
  "BilaplacianEigs_ShiftedInverseEigenvalue":
    ("bilaplacian_convergence_eigenvalue_shifted_inverse", "eigenvalue_approx_SI", (1, 1, 2)),
  "Diffusion_MassApproxEigenvalue":
    ("diffusion_convergence_eigenvalue_approx", "eigenvalue_approx_MA", (1, 1, 2)),
  "Bilaplacian_MassApproxEigenvalue":
    ("bilaplacian_convergence_eigenvalue_approx", "eigenvalue_approx_MA", (1, 1, 2)),
}

## Uniform bilaplacian solver, for which no reproducible exists.
def bilaplacian_uniform():
  import numpy as np
  import scipy.sparse.linalg as sp_lin_alg
  import HyperHDG
  const                     = HyperHDG.config()
  const.global_loop         = "Elliptic"
  const.local_solver        = "BilaplacianUniform < 1, 1, 2 >"
  const.topology            = "Cubic< 1, 2 >"
  const.geometry            = "UnitCube< 1, 2 >"
  const.node_descriptor     = "Cubic< 1, 2 >"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
  hdg_wrapper = HyperHDG.include(const)([3, 3])
  vector_dirichlet = hdg_wrapper.zero_vector()
  vector_dirichlet[0] = 1.
  hdg_wrapper.read_dirichlet_indices(np.array([0, 1, len(vector_dirichlet) - 2,
                                               len(vector_dirichlet) - 1]))
  rhs = [-i for i in hdg_wrapper.trace_to_flux(vector_dirichlet)]
  n = hdg_wrapper.size_of_system()
  A = sp_lin_alg.LinearOperator((n, n), matvec=hdg_wrapper.trace_to_flux)
  solution, info = sp_lin_alg.gmres(A, rhs, tol=1e-10, maxiter=1000)
  assert info == 0 and np.all(np.isfinite(solution)), "Linear solver did not converge!"

if __name__ == "__main__":
  assert len(sys.argv) == 2, "Usage: python3 run.py <case name>"
  case = sys.argv[1]
  if case == "BilaplacianUniform_Elliptic":
    bilaplacian_uniform()
  else:
    module, function, args = CASES[case]
    getattr(importlib.import_module(module), function)(*args)
  print("Local solver test", case, "succeeded.")
