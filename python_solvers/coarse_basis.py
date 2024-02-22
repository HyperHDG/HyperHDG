import numpy as np
import scipy.sparse as sp


def _coarse_basis(points, n_elem_1d, h=None, epsilon=1e-10):
  # B contains (N+1)^2 columns and length(p) rows.
  # Each column represents a bilinear basis function on the mesh of N*N squares of side length H=1/N
  if h is None:       h = 1. / n_elem_1d

  m, n = (n_elem_1d+1)**2, len(points)  # Dimensions of the resulting matrix
  node_vec, index_vec, value_vec = [], [], []

  # loop over elements to find nodes
  for i in range(n_elem_1d):
    for j in range(n_elem_1d):
      nodes = np.where( (points[:,0] + epsilon > h * i) & (points[:,0] + epsilon < h * i + h) &
                        (points[:,1] + epsilon > h * j) & (points[:,1] + epsilon < h * j + h) )[0]
      index = [ j*(n_elem_1d+1)+i,       j*(n_elem_1d+1)+i+1, 
                (j+1)*(n_elem_1d+1)+i+1, (j+1)*(n_elem_1d+1)+i ]
      x, y  = points[nodes,0] / h - i, points[nodes,1] / h - j

      node_vec.extend( [nodes, nodes, nodes, nodes] )
      index_vec.extend( [index[0]*np.ones(len(nodes)), index[1]*np.ones(len(nodes)),
                         index[2]*np.ones(len(nodes)), index[3]*np.ones(len(nodes))] )
      value_vec.extend( [(1-x)*(1-y), x*(1-y), x*y, (1-x)*y] )

  node_vec = np.hstack(node_vec)
  index_vec = np.hstack(index_vec)
  value_vec = np.hstack(value_vec)

  return sp.csr_matrix((value_vec, (node_vec, index_vec)), shape=(n, m))


class precond_data:
  def __init__( self, points, n_elem_1d, h=None, epsilon=1e-10 ):
    coarse_basis = _coarse_basis(points, n_elem_1d, h, epsilon)

    n_elem_1d = int(np.sqrt(coarse_basis.shape[1]) - 1)
    int_nodes = np.concatenate([ j*(n_elem_1d+1) + np.arange(2, n_elem_1d+1) - 1 
                                 for j in range(1, n_elem_1d) ])  # find interior nodes
    coarse_basis_int = coarse_basis[:, int_nodes]
    bnd_nodes = np.where(np.sum(coarse_basis_int, axis=1) < epsilon)[0]  # coarse basis in V

    # Get the row indices and data of the column to be set to zero
    for index in bnd_nodes:
      start_idx, end_idx = coarse_basis.indptr[index], coarse_basis.indptr[index + 1]
      coarse_basis.data[start_idx:end_idx] = 0.

    self.coarse_basis     = coarse_basis
    self.coarse_basis_int = coarse_basis_int


def precond(lhs_mat, precond_data, rhs_vec, epsilon=1e-8):
  coarse_basis = precond_data.coarse_basis
  coarse_basis_int = precond_data.coarse_basis_int

  result_mat = coarse_basis_int.dot(
               sp.linalg.spsolve(coarse_basis_int.T.dot(lhs_mat.dot(coarse_basis_int)),
               coarse_basis_int.T.dot(rhs_vec)) )  # coarse solve

  for k in range(coarse_basis.shape[1]):
    nj = np.nonzero(coarse_basis.getcol(k) > epsilon)[0]
    Ij = np.zeros((coarse_basis.shape[0], len(nj)))
    Ij[nj, np.arange(len(nj))] = 1
    result_mat += Ij.dot(sp.linalg.spsolve(lhs_mat[nj, :][:, nj], Ij.T.dot(rhs_vec)))  # local solve
  
  return result_mat
