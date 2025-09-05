#!/usr/bin/env python3

import numpy as np
import scipy.sparse as sp
import sys
import datetime

# TODO: optimize the coarse_basis creation / swap it out for algebraic construction

def _coarse_basis_2d(points, n_elem_1d, epsilon=1e-10):
  # B contains (N+1)^2 columns and length(p) rows.
  # Each column represents a bilinear basis function on the mesh of N*N squares of side length H=1/N

  min_coord = np.min(points, axis=0)
  h = np.divide( (np.max(points, axis=0) - min_coord)[:len(n_elem_1d)], n_elem_1d )

  m, n = (n_elem_1d[0]+1) * (n_elem_1d[1]+1), len(points)  # Dimensions of the resulting matrix
  node_vec, index_vec, value_vec = [], [], []

  # loop over elements to find nodes
  for i in range(n_elem_1d[0]):
    for j in range(n_elem_1d[1]):
      nodes = np.where( (points[:,0] - min_coord[0] + epsilon > h[0] * i) &
                        (points[:,0] - min_coord[0] + epsilon < h[0] * (i + 1)) &
                        (points[:,1] - min_coord[1] + epsilon > h[1] * j) &
                        (points[:,1] - min_coord[1] + epsilon < h[1] * (j + 1)) )[0]
      index = [ j*(n_elem_1d[0]+1)+i,       j*(n_elem_1d[0]+1)+i+1,
                (j+1)*(n_elem_1d[0]+1)+i+1, (j+1)*(n_elem_1d[0]+1)+i ]
      x = (points[nodes,0] - min_coord[0]) / h[0] - i
      y = (points[nodes,1] - min_coord[1]) / h[1] - j

      node_vec.extend( [nodes, nodes, nodes, nodes] )
      index_vec.extend( [index[0]*np.ones(len(nodes)), index[1]*np.ones(len(nodes)),
                         index[2]*np.ones(len(nodes)), index[3]*np.ones(len(nodes))] )
      value_vec.extend( [(1-x)*(1-y), x*(1-y), x*y, (1-x)*y] )

  node_vec, index_vec, value_vec = np.hstack(node_vec), np.hstack(index_vec), np.hstack(value_vec)
  return sp.csr_matrix((value_vec, (node_vec, index_vec)), shape=(n, m))


def _coarse_basis_3d(points, n_elem_1d, epsilon=1e-10):
  # B contains (N+1)^3 columns and length(p) rows.
  # Each column represents a bilinear basis function on the mesh of N*N squares of side length H=1/N

  min_coord = np.min(points, axis=0)
  h = np.divide( (np.max(points, axis=0) - min_coord)[:len(n_elem_1d)], n_elem_1d )

  m, n = (n_elem_1d[0]+1) * (n_elem_1d[1]+1) * (n_elem_1d[2]+1), len(points)
  node_vec, index_vec, value_vec = [], [], []

  # loop over elements to find nodes
  for i in range(n_elem_1d[0]):
    for j in range(n_elem_1d[1]):
      for k in range(n_elem_1d[2]):
        nodes = np.where( (points[:,0] - min_coord[0] + epsilon > h[0] * i) &
                          (points[:,0] - min_coord[0] + epsilon < h[0] * (i + 1)) &
                          (points[:,1] - min_coord[1] + epsilon > h[1] * j) &
                          (points[:,1] - min_coord[1] + epsilon < h[1] * (j + 1)) &
                          (points[:,2] - min_coord[2] + epsilon > h[2] * k) &
                          (points[:,2] - min_coord[2] + epsilon < h[2] * (k + 1)) )[0]
        index = [ (k * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i,
                  (k * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i+1,
                  (k * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i+1,
                  (k * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i,
                  ((k+1) * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i,
                  ((k+1) * (n_elem_1d[1]+1) + j) * (n_elem_1d[0]+1) + i+1,
                  ((k+1) * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i+1,
                  ((k+1) * (n_elem_1d[1]+1) + j+1) * (n_elem_1d[0]+1) + i ]
        x = (points[nodes,0] - min_coord[0]) / h[0] - i
        y = (points[nodes,1] - min_coord[1]) / h[1] - j
        z = (points[nodes,2] - min_coord[2]) / h[2] - k

        node_vec.extend( [nodes, nodes, nodes, nodes, nodes, nodes, nodes, nodes] )
        index_vec.extend( [index[0]*np.ones(len(nodes)), index[1]*np.ones(len(nodes)),
                           index[2]*np.ones(len(nodes)), index[3]*np.ones(len(nodes)),
                           index[4]*np.ones(len(nodes)), index[5]*np.ones(len(nodes)),
                           index[6]*np.ones(len(nodes)), index[7]*np.ones(len(nodes))] )
        value_vec.extend( [ (1-z)*(1-x)*(1-y), (1-z)*x*(1-y), (1-z)*x*y, (1-z)*(1-x)*y,
                            z*(1-x)*(1-y), z*x*(1-y), z*x*y, z*(1-x)*y ] )

  node_vec, index_vec, value_vec = np.hstack(node_vec), np.hstack(index_vec), np.hstack(value_vec)
  return sp.csr_matrix((value_vec, (node_vec, index_vec)), shape=(n, m))


def coarse_basis_pu(points, domains):
  # nodal values of all the basis functions combined
  nodal = np.zeros(len(points))
  for domain in domains:
    nodal[domain] += 1.
  mask = nodal != 0
  nodal[mask] = 1/nodal[mask]

  # now split the nodal values to each of the basis functions
  # in the overlap, the nodal values are shared
  nnzs = np.zeros(len(domains.all_domains))
  for k in range(len(domains)):
    start = domains.ioffsets[k]
    end = domains.ioffsets[k+1]
    nnzs[start:end] = nodal[domains.all_domains[start:end]]

  return sp.csc_matrix((nnzs, domains.all_domains, domains.ioffsets))


class JPrecond:
  def __init__( self, lhs_mat, points, n_elem_1d, epsilon=1e-10, repeat=1, domains=None, coarse_space=None):
    # save lhs_mat
    self.lhs_mat = lhs_mat
    self.coarse_space = coarse_space

    coarse_basis, int_nodes = [], []
    if coarse_space.lower() == "pu":
      coarse_basis_int = coarse_basis_pu(points, domains)
    elif coarse_space.lower() == "q1":
      if   len(n_elem_1d) == 2:
        coarse_basis = _coarse_basis_2d(points, n_elem_1d,  epsilon)
        int_nodes    = np.concatenate([ j*(n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
                                        for j in range(1, n_elem_1d[1]) ])  # find interior nodes
        coarse_basis_int = coarse_basis[:, int_nodes]
      elif len(n_elem_1d) == 3:
        coarse_basis = _coarse_basis_3d(points, n_elem_1d, epsilon)
        int_nodes    = np.concatenate([
            (k * (n_elem_1d[1]+1)+j) * (n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
            for j in range(1, n_elem_1d[1]) for k in range(1, n_elem_1d[2]) ])  # find interior nodes
        coarse_basis_int = coarse_basis[:, int_nodes]
      else:
          raise RuntimeError(f"expected len(n_elem_1d)=2,3 got {len(n_elem_1d)}")
    else:
      raise RuntimeError(f"unexpected coarse space '{coarse_space}'")

    if repeat > 1:
      coarse_basis_int = sp.kron(coarse_basis_int, np.eye(repeat))

    self.coarse_basis_int = sp.csc_matrix(coarse_basis_int)

    # NOTE: matrix-free: self.coarse_basis_int is of size (n,m) where m << n,
    #       hence, product is still faster than assembly of the full system
    # precond lhs_mat and explicitly format in csc
    precond_lhs = sp.csc_matrix(
      self.coarse_basis_int.T @ self.lhs_mat @ self.coarse_basis_int
    )
    # precompute splu of precond_lhs
    start = datetime.datetime.now()
    self.splu_precond_lhs = sp.linalg.splu(precond_lhs)
    self.init_clu_time = (datetime.datetime.now() - start).total_seconds()

    # NOTE: matrix-free: the submatrices essentially form an overlapping block-diagonal submatrix
    #       hence, could optimize to only construct that
    # precompute splu
    ioffsets_r = repeat*domains.ioffsets
    all_domains_r = repeat*np.repeat(domains.all_domains, repeat) + np.tile(np.arange(repeat), len(domains.all_domains))
    self.domains = Domains(ioffsets_r, all_domains_r)

    start = datetime.datetime.now()
    submats = [self.lhs_mat[nj, :][:, nj] for nj in self.domains]
    self.init_submat_time = (datetime.datetime.now() - start).total_seconds()

    start = datetime.datetime.now()
    self.splu = [sp.linalg.splu(submat) for submat in submats]
    self.init_lu_time = (datetime.datetime.now() - start).total_seconds()

    self.total_solve_time = 0
    self.total_csolve_time = 0

  def matmul(self, rhs_vec, epsilon=1e-14):
    start = datetime.datetime.now()
    precond_rhs = self.coarse_basis_int.T @ rhs_vec
    result_vec = self.coarse_basis_int @ self.splu_precond_lhs.solve(precond_rhs)
    self.total_csolve_time += (datetime.datetime.now() - start).total_seconds()

    start = datetime.datetime.now()
    for nj, splu in zip(self.domains, self.splu):
      result_vec[nj] += splu.solve(rhs_vec[nj])
    self.total_solve_time += (datetime.datetime.now() - start).total_seconds()

    return result_vec

class Domains:
  def __init__(self, ioffsets, all_domains):
    self.ioffsets = ioffsets
    self.all_domains = all_domains

  def __getitem__(self, dom_index):
    if dom_index+1 >= len(self.ioffsets):
      raise IndexError()

    start = self.ioffsets[dom_index]
    end = self.ioffsets[dom_index+1]
    return self.all_domains[start:end]

  def __len__(self):
    return len(self.ioffsets)-1
