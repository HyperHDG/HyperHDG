# --------------------------------------------------------------------------------------------------
# The data necessary to run this file can be obtained using the HyperHDG.fiber_network.make_geo
# function on the files of Hauck, M., & Rupp, A. (2024). Fiber network models of paper. Zenodo.
# https://doi.org/10.5281/zenodo.12751486
# --------------------------------------------------------------------------------------------------

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import os, sys, argparse, logging, datetime, subprocess

####### PRECONDITIONER CODE

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


class JPrecond:
  def __init__( self, lhs_mat, points, n_elem_1d, epsilon=1e-10, repeat=1 ):
    # save lhs_mat
    self.lhs_mat = lhs_mat

    coarse_basis, int_nodes = [], []
    if   len(n_elem_1d) == 2:
      coarse_basis = _coarse_basis_2d(points, n_elem_1d,  epsilon)
      int_nodes    = np.concatenate([ j*(n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
                                      for j in range(1, n_elem_1d[1]) ])  # find interior nodes
    elif len(n_elem_1d) == 3:
      coarse_basis = _coarse_basis_3d(points, n_elem_1d, epsilon)
      int_nodes    = np.concatenate([
        (k * (n_elem_1d[1]+1)+j) * (n_elem_1d[0]+1) + np.arange(1, n_elem_1d[0])
        for j in range(1, n_elem_1d[1]) for k in range(1, n_elem_1d[2]) ])  # find interior nodes

    coarse_basis_int = coarse_basis[:, int_nodes]
    bnd_nodes = np.where(np.sum(coarse_basis_int, axis=1) < epsilon)[0]  # coarse basis in V

    # Get the row indices and data of the column to be set to zero
    for index in bnd_nodes:
      start_idx, end_idx = coarse_basis.indptr[index], coarse_basis.indptr[index + 1]
      coarse_basis.data[start_idx:end_idx] = 0.

    if repeat > 1:
      coarse_basis = sp.kron(coarse_basis, np.eye(repeat))
      coarse_basis_int = sp.kron(coarse_basis_int, np.eye(repeat))

    self.coarse_basis     = sp.csc_matrix(coarse_basis)
    self.coarse_basis_int = sp.csc_matrix(coarse_basis_int)

    # precompute splu of precond_lhs
    self.splu_precond_lhs = sp.linalg.splu(self.coarse_basis_int.T @ self.lhs_mat @ self.coarse_basis_int)

    # precompute splu
    self.splu = [None] * self.coarse_basis.shape[1]
    for k in range(self.coarse_basis.shape[1]):
      col_k = self.coarse_basis.getcol(k)
      nj = col_k.indices[col_k.data > epsilon]
      if nj.size == 0:
        continue
      self.splu[k] = sp.linalg.splu(self.lhs_mat[nj, :][:, nj])

  def matmul(self, rhs_vec, epsilon=1e-14):
    precond_rhs = self.coarse_basis_int.T @ rhs_vec
    result_vec = self.coarse_basis_int @ self.splu_precond_lhs.solve(precond_rhs)

    for k in range(self.coarse_basis.shape[1]):
      col_k = self.coarse_basis.getcol(k)
      nj = col_k.indices[col_k.data > epsilon]
      if nj.size == 0:
        continue
      result_vec[nj] += self.splu[k].solve(rhs_vec[nj])
    return result_vec


######### SETUP argument parsing & logging

time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
git_hash =  subprocess.run(["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()

logger = logging.getLogger("fiber_network_elastic")
logger.setLevel(logging.INFO)
log_formatter  = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

os.system("mkdir -p logs")
fhandler = logging.FileHandler(f"logs/fiber_network_elastic.{time_stamp}.log", mode="w")
fhandler.setFormatter(log_formatter)
fhandler.setLevel(logging.INFO)

shandler = logging.StreamHandler(sys.stdout)
shandler.setFormatter(log_formatter)
shandler.setLevel(logging.INFO)

logger.addHandler(fhandler)
logger.addHandler(shandler)

parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("network")
parser.add_argument("-t", "--rtol",
  help="relative tolerance when to stop the CG iterator",
  type=float, default=1e-10
)

args = parser.parse_args()
domain = args.network

######## MAIN code

logger.info("started with")
logger.info(f"  git_hash={git_hash}")
logger.info(f"  network={args.network}")
logger.info(f"  rtol={args.rtol}")

os.system("mkdir -p output")

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG
  
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.local_solver    = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP( domain + ".geo" )
domain_points = np.loadtxt(domain + "_points.txt")

rhs = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

logger.info("assembling  A...")

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

logger.info("assembling  B...")

precond = JPrecond(A, domain_points, [2**3, 2**3], repeat=6)
B = sp.linalg.LinearOperator(
  (system_size,system_size),
  matvec=precond.matmul
)

iters = 0
start = datetime.datetime.now()

def log_iter(x):
  global iters, start
  iters += 1
  relErr = np.linalg.norm(A @ x - rhs) / np.linalg.norm(rhs)
  bilin = .5 * x.dot(A @ x) - x.dot(rhs)
  duration = datetime.datetime.now() - start
  start = datetime.datetime.now()
  logger.info(f"{iters:>5} {relErr:>13.6e} {bilin:>13.6e} {duration}")

logger.info("starting cg")

vectorSolution, num_iter = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=log_iter, M=B)

if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

error = HDG_wrapper.errors(vectorSolution)[0]
logger.info(f"HDG_wrapper error={error:>.6e}")

output_name = domain + "_timo"
HDG_wrapper.plot_option("fileName", output_name)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution)

logger.info(f"solution written to 'output/{output_name}'")
