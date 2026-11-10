# Check GlobalLoop::GeneralizedEigenvalue with LocalSolver::DiffusionEigsWeighted.
#
# - Without weights and hypernode mass, the Schur complement of the stiffness matrix with respect to
#   the bulk unknowns is (minus) the condensed matrix of LocalSolver::Diffusion.
# - The hyperedge properties (d, kappa, rho, c, mu) of a .geo file scale the local matrices as
#   intended, and the hypernode mass acts on the faces of node type 3 only.
# - The eigenvalues of the Dirichlet Laplacian on the unit square converge with order 2p + 1.
# - The evaluation helpers (bulk_values, face_values, ...) integrate discrete functions exactly.

import os, sys, tempfile
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../python")
import hyperhdg

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/topology/file.hxx>

#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/global_loop/generalized_eigenvalue.hxx>
#include <HyperHDG/local_solver/diffusion_eigs_weighted_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>

template <unsigned int space_dimT, typename param_float_t = double>
struct BoundaryParameters
{{
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{{1, 2, 3, 4, 5, 6, 7, 8}};
}};

using Cubic = GlobalLoop::GeneralizedEigenvalue<
  Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>, NodeDescriptor::Cubic<2, 2>,
  LocalSolver::DiffusionEigsWeighted<2, {degree}, {quadrature}, BoundaryParameters> >;
using CubicNoDirichlet = GlobalLoop::GeneralizedEigenvalue<
  Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>, NodeDescriptor::Cubic<2, 2>,
  LocalSolver::DiffusionEigsWeighted<2, {degree}, {quadrature}, LocalSolver::DiffusionParametersDefault> >;
using Condensed = GlobalLoop::Elliptic<
  Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>, NodeDescriptor::Cubic<2, 2>,
  LocalSolver::Diffusion<2, {degree}, {quadrature}, LocalSolver::DiffusionParametersDefault> >;
using File = GlobalLoop::GeneralizedEigenvalue<
  Topology::File<2, 3>, Geometry::File<2, 3>, NodeDescriptor::File<2, 3>,
  LocalSolver::DiffusionEigsWeighted<2, {degree}, {quadrature}> >;

NB_MODULE(hdg_test_generalized_eigenvalue, m) {{
  HyperHDG::bind_python<Cubic, std::vector<unsigned int> >(m, "Cubic");
  HyperHDG::bind_python<CubicNoDirichlet, std::vector<unsigned int> >(m, "CubicNoDirichlet");
  HyperHDG::bind_python<Condensed, std::vector<unsigned int> >(m, "Condensed");
  HyperHDG::bind_python<File>(m, "File");
}}
"""

def matrices(hdg):
  n = hdg.size_of_system()
  rows, cols, vals = hdg.stiffness_mat()
  stiff = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
  rows, cols, vals = hdg.mass_mat()
  mass = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
  return stiff, mass

def eigenvalues(hdg, k):
  stiff, mass = matrices(hdg)
  keep = np.setdiff1d(np.arange(stiff.shape[0]), np.asarray(hdg.dirichlet_indices(), dtype=int))
  stiff, mass = stiff[keep][:, keep].tocsc(), mass[keep][:, keep].tocsc()
  rank = int(np.count_nonzero(abs(mass).sum(axis=1).A1))
  vals = sp_lin_alg.eigsh(stiff, k=k, M=mass, sigma=0., which="LM", ncv=min(20, rank),
                          return_eigenvectors=False)
  return np.sort(vals)

def write_plate(filename, props):
  # one unit square in the plane z = 0 of R^3: faces x = 0 (mass), x = 1 (Dirichlet), y = 0, 1
  with open(filename, "w") as file:
    file.write("Space_Dim = 3;\nHyperEdge_Dim = 2;\nN_Points = 4;\nN_HyperNodes = 4;\n"
               "N_HyperEdges = 1;\nPOINTS:\n0 0 0\n1 0 0\n0 1 0\n1 1 0\n"
               "HYPERNODES_OF_HYPEREDGES:\n0 1 2 3\nTYPES_OF_HYPERFACES:\n3 1 2 2\n"
               "POINTS_OF_HYPEREDGES:\n0 1 2 3\n")
    if props is not None:
      file.write("HYPEREDGE_PROPERTIES: %d\n%s\n" % (len(props), " ".join(map(repr, props))))

success = True
degree = 2
module = hyperhdg.load(MODULE.format(degree=degree, quadrature=2 * degree + 2))

# Schur complement = - condensed matrix of LocalSolver::Diffusion (tau = 0.7).
hdg, condensed = module.CubicNoDirichlet([3, 3], 0.7), module.Condensed([3, 3], 0.7)
stiff = matrices(hdg)[0].toarray()
ns = hdg.n_skeleton_dofs()
schur = stiff[:ns, :ns] - stiff[:ns, ns:] @ np.linalg.solve(stiff[ns:, ns:], stiff[ns:, :ns])
rows, cols, vals = condensed.trace_to_flux_mat()
reference = sp.coo_matrix((vals, (rows, cols)), shape=(ns, ns)).toarray()
deviation = np.abs(schur + reference).max() / np.abs(reference).max()
ok = deviation < 1e-12 and np.abs(stiff - stiff.T).max() < 1e-12 * np.abs(stiff).max()
success &= ok
print("Schur complement vs. condensed Diffusion matrix: deviation %.1e %s"
      % (deviation, "OK" if ok else "FAILED"))

# Hyperedge properties (d, kappa, rho, c, mu) and hypernode mass.
with tempfile.TemporaryDirectory() as directory:
  plain, weighted = os.path.join(directory, "plain.geo"), os.path.join(directory, "weighted.geo")
  write_plate(plain, None)
  write_plate(weighted, [2., 1., 3., 0., 5.])
  a0, m0 = matrices(module.File(plain, 1.))
  hdg = module.File(weighted, 1.)
  a1, m1 = matrices(hdg)
  nb, ns, nf = hdg.n_loc_bulk_dofs(), hdg.n_skeleton_dofs(), degree + 1
  mass_face = m1[:ns, :ns].toarray()
  ok = abs(a1 - 2. * a0).max() < 1e-12 and abs(m1[ns:, ns:] - 6. * m0[ns:, ns:]).max() < 1e-12 \
       and abs(m0[:ns, :ns]).max() == 0. \
       and np.allclose(mass_face[:nf, :nf], 5. * np.eye(nf)) and abs(mass_face[nf:, :]).max() == 0. \
       and np.allclose(hdg.coefficients(), [2., 1., 3., 0., 5.])
  # c = 1 adds d * (u, v) to the stiffness matrix
  write_plate(weighted, [2., 1., 3., 1., 5.])
  a2 = matrices(module.File(weighted, 1.))[0]
  ok = ok and abs((a2 - a1)[ns:, ns:] - 2. * m0[ns:, ns:]).max() < 1e-12 \
       and abs((a2 - a1)[:ns, :]).max() < 1e-12
  # quadrature with the evaluation helpers: int_E u = (u, 1)_E and int_N u_hat = <u_hat, 1>_N
  x = np.random.default_rng(1).standard_normal(hdg.size_of_system())
  points, weights = np.polynomial.legendre.leggauss(degree + 2)
  points, weights = 0.5 * (points + 1.), 0.5 * weights
  bulk = hdg.bulk_values(list(x), list(points)) @ np.outer(weights, weights).ravel()
  face = hdg.face_values(list(x), list(points)).reshape(4, -1)[0] @ weights
  ok = ok and abs(bulk - x[ns]) < 1e-12 and abs(face - x[0]) < 1e-12 \
       and np.allclose(hdg.face_node_types(), [3, 1, 2, 2]) \
       and np.allclose(hdg.bulk_points([0.5]), [0.5, 0.5, 0.]) \
       and np.allclose(hdg.face_points([0.25]).reshape(4, 3)[0], [0., 0.25, 0.])
success &= ok
print("hyperedge properties, hypernode mass, evaluation:", "OK" if ok else "FAILED")

# Dirichlet eigenvalues of the unit square: 2 pi^2, 5 pi^2, 5 pi^2 with order 2p + 1.
exact = np.pi ** 2 * np.array([2., 5., 5.])
errors = [np.abs(eigenvalues(module.Cubic([n, n], 1.), 3) - exact).max() for n in [8, 16]]
rate = np.log2(errors[0] / errors[1])
ok = rate > 2 * degree + 0.5
success &= ok
print("eigenvalue errors", errors, "order %.2f" % rate, "OK" if ok else "FAILED")

assert success, "Some checks of GlobalLoop::GeneralizedEigenvalue failed!"
