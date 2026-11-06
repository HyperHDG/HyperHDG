"""Shared helpers of the thin-plate eigenvalue experiments.

- compile and load the HyperHDG module (GlobalLoop::GeneralizedEigenvalue with
  LocalSolver::DiffusionEigsWeighted) through python/hyperhdg.py,
- write hypergraphs of two-dimensional hyperedges to .geo files,
- assemble the stiffness and mass matrices and compute eigenpairs with shift and invert,
- Gauss quadrature, the gap between subspaces, and experimental orders of convergence.
"""

import math
import os
import sys

import numpy as np
import scipy.linalg
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
# Generated domain files (output/ is ignored by git).
OUTPUT = os.path.join(ROOT, "output", "thin_plate_eigenvalues")
sys.path.insert(0, os.path.join(ROOT, "python"))
import hyperhdg  # noqa: E402  (python/hyperhdg.py of this HyperHDG tree)

# Node types of the faces in the .geo files, cf. LocalSolver::DiffusionEigsWeightedParametersDefault.
INTERIOR, DIRICHLET, NATURAL, MASS = 0, 1, 2, 3

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <HyperHDG/global_loop/generalized_eigenvalue.hxx>
#include <HyperHDG/local_solver/diffusion_eigs_weighted_ldgh.hxx>

using HDG = GlobalLoop::GeneralizedEigenvalue<
  Topology::File<2, {space_dim}>,
  Geometry::File<2, {space_dim}>,
  NodeDescriptor::File<2, {space_dim}>,
  LocalSolver::DiffusionEigsWeighted<2, {degree}, {quadrature}>,
  std::vector<double> >;

NB_MODULE(thin_plate_eigenvalues, m) {{ HyperHDG::bind_python<HDG>(m, "{name}"); }}
"""


def problem_class(degree, space_dim):
  """Compile (once, cached afterwards) and return the class for degree and space dimension."""
  name = f"GeneralizedEigenvalue_S{space_dim}_P{degree}"
  code = MODULE.format(space_dim=space_dim, degree=degree, quadrature=2 * degree + 2, name=name)
  cache = os.environ.get("HYPERHDG_CACHE_DIR", os.path.join(ROOT, ".hyperhdg-cache"))
  return hyperhdg.load(code, cache_dir=cache)[name]


def write_geo(filename, space_dim, points, edges, comment=""):
  """Write a hypergraph of two-dimensional hyperedges (parallelograms) to a .geo file.

  points: array of shape (n_points, space_dim).
  edges:  list of dicts with the keys
          'points' -- four point indices in lexicographic order (p0, p0 + a, p0 + b, p0 + a + b),
          'nodes'  -- four hypernode indices of the faces (a = 0, a = 1, b = 0, b = 1),
          'types'  -- four node types of these faces,
          'props'  -- hyperedge properties (d, kappa, rho, c, mu), cf. DiffusionEigsWeighted.
  """
  n_nodes = 1 + max(max(edge["nodes"]) for edge in edges)
  n_props = len(edges[0]["props"])
  with open(filename, "w") as file:
    if comment:
      file.write("".join(f"# {line}\n" for line in comment.splitlines()))
    file.write(f"\nSpace_Dim     = {space_dim};\nHyperEdge_Dim = 2;\n\n")
    file.write(f"N_Points      = {len(points)};\nN_HyperNodes  = {n_nodes};\n")
    file.write(f"N_HyperEdges  = {len(edges)};\n\nPOINTS:\n")
    for point in points:
      file.write(" ".join(repr(float(coord)) for coord in point) + "\n")
    file.write("\nHYPERNODES_OF_HYPEREDGES:\n")
    file.write("".join(" ".join(map(str, edge["nodes"])) + "\n" for edge in edges))
    file.write("\nTYPES_OF_HYPERFACES:\n")
    file.write("".join(" ".join(map(str, edge["types"])) + "\n" for edge in edges))
    file.write("\nPOINTS_OF_HYPEREDGES:\n")
    file.write("".join(" ".join(map(str, edge["points"])) + "\n" for edge in edges))
    file.write(f"\nHYPEREDGE_PROPERTIES: {n_props}\n")
    file.write("".join(" ".join(repr(float(v)) for v in edge["props"]) + "\n" for edge in edges))


class Discretization:
  """Assembled generalized eigenvalue problem A x = lambda M x without the Dirichlet unknowns."""

  def __init__(self, hdg):
    self.hdg = hdg
    self.size = hdg.size_of_system()
    rows, cols, vals = hdg.stiffness_mat()
    stiff = sp.csr_matrix((vals, (rows, cols)), shape=(self.size, self.size))
    rows, cols, vals = hdg.mass_mat()
    mass = sp.csr_matrix((vals, (rows, cols)), shape=(self.size, self.size))
    dirichlet = np.asarray(hdg.dirichlet_indices(), dtype=np.int64)
    self.keep = np.setdiff1d(np.arange(self.size), dirichlet)
    self.stiff = stiff[self.keep][:, self.keep].tocsc()
    self.mass = mass[self.keep][:, self.keep].tocsc()
    self.mass_full = mass

  def eigenpairs(self, n_eigs, sigma=0., refinement_steps=2, extra=4, tol=1e-12):
    """The n_eigs eigenpairs closest to sigma, sorted, with M-orthonormal full-length vectors.

    The problem is equilibrated by the symmetric diagonal scaling with diag(A)^(-1/2), which does
    not change the eigenvalues, and the shifted systems are solved by a sparse LU factorization
    with iterative refinement. Both matter for strongly graded and anisotropic meshes, whose
    matrices have entries of very different magnitude.

    Lanczos (ARPACK in shift-and-invert mode) computes n_eigs + extra eigenpairs, since the last
    ones of the requested set are unreliable if they split a cluster of eigenvalues, and starts
    from a fixed pseudo-random vector, so that runs are reproducible. Since M is singular, the
    computed vectors are purified by one step of inverse iteration, which removes spurious
    components in the kernel of M, and the eigenpairs are the Ritz pairs of the purified vectors.
    """
    scaling = sp.diags(1. / np.sqrt(abs(self.stiff.diagonal())))
    stiff = (scaling @ self.stiff @ scaling).tocsc()
    mass = (scaling @ self.mass @ scaling).tocsc()
    shifted = (stiff - sigma * mass).tocsc()
    lu = spla.splu(shifted)

    def solve(rhs):
      x = lu.solve(rhs)
      for _ in range(refinement_steps):
        x += lu.solve(rhs - shifted @ x)
      return x

    n = stiff.shape[0]
    op_inv = spla.LinearOperator((n, n), matvec=solve, dtype=float)
    rank = int(np.count_nonzero(abs(mass).sum(axis=1).A1))
    k = min(n_eigs + extra, rank - 1)
    ncv = min(max(2 * k + 1, 20), rank)
    start = np.random.default_rng(0).standard_normal(n)
    vals, vecs = spla.eigsh(stiff, k=k, M=mass, sigma=sigma, which="LM", ncv=ncv, tol=tol,
                            OPinv=op_inv, v0=start)
    vecs = vecs[:, np.argsort(vals)]
    vecs = np.column_stack([solve(mass @ vecs[:, i]) for i in range(k)])  # purification
    ritz_vals, coeffs = scipy.linalg.eigh(vecs.T @ (stiff @ vecs), vecs.T @ (mass @ vecs))
    chosen = np.sort(np.argsort(abs(ritz_vals - sigma))[:n_eigs])  # closest to sigma, ascending
    vals, vecs = ritz_vals[chosen], scaling @ (vecs @ coeffs[:, chosen])
    full = np.zeros((self.size, n_eigs))
    full[self.keep] = vecs  # M-orthonormal, since the Ritz vectors are
    self.residuals = [np.linalg.norm(self.stiff @ full[self.keep, i] - vals[i] *
                                     (self.mass @ full[self.keep, i])) /
                      np.linalg.norm(self.stiff @ full[self.keep, i]) for i in range(n_eigs)]
    return vals, full

  def n_unknowns(self):
    return len(self.keep)


def gauss(n):
  """Gauss-Legendre points and weights on [0, 1]."""
  points, weights = np.polynomial.legendre.leggauss(n)
  return 0.5 * (points + 1.), 0.5 * weights


def gap(gram_u, gram_v, cross):
  """Gap between span(u_i) and span(v_j) from the Gram matrices and cross[i, j] = <u_i, v_j>.

  For subspaces of equal dimension, both directed gaps coincide with the sine of the largest
  principal angle. If the dimensions differ, the larger directed gap (equal to one) is returned.
  """
  inv_u = np.linalg.inv(np.linalg.cholesky(gram_u))
  inv_v = np.linalg.inv(np.linalg.cholesky(gram_v))
  cosines = np.linalg.svd(inv_u @ cross @ inv_v.T, compute_uv=False)
  if gram_u.shape[0] != gram_v.shape[0]:
    return 1.
  return math.sqrt(max(0., 1. - min(1., cosines.min()) ** 2))


def eoc(errors, sizes):
  """Experimental orders of convergence log(e_{k-1}/e_k) / log(h_{k-1}/h_k), None for the first."""
  orders = [None]
  for k in range(1, len(errors)):
    if errors[k] > 0. and errors[k - 1] > 0.:
      orders.append(math.log(errors[k - 1] / errors[k]) / math.log(sizes[k - 1] / sizes[k]))
    else:
      orders.append(None)
  return orders
