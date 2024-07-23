from __future__ import print_function

import numpy as np
import scipy as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_approx():
  def __init__(self, hdg_wrapper):
    self.hdg_wrapper       = hdg_wrapper
    dirichlet_indices = hdg_wrapper.dirichlet_nodes()
    self.dirichlet_indices = [index * 6 + comp for index in dirichlet_indices for comp in range(6)]
    self.index_vector      = [-1] * self.short_vector_size()
    self.vector_index      = [-1] * self.long_vector_size()
    n_indices = 0
    for i in range(len(self.index_vector)):
      while n_indices < len(self.dirichlet_indices) and \
            i + n_indices == self.dirichlet_indices[n_indices]:
        n_indices = n_indices + 1
      self.index_vector[i] = i + n_indices
      self.vector_index[i + n_indices] = i
    assert -1 not in self.index_vector
    self.index_dict = {value: idx for idx, value in enumerate(self.index_vector)}
  def long_vector_size(self):
    return self.hdg_wrapper.size_of_system()
  def short_vector_size(self):
    return self.hdg_wrapper.size_of_system() - len(self.dirichlet_indices)
  def long_vector(self, vector):
    return [vector[x] if x > -1 else 0. for x in self.vector_index]
  def short_vector(self, vector):
    return [vector[x] for x in self.index_vector]
  def stiff_mat(self, eig = 0.):
    system_size = self.short_vector_size()
    col_ind, row_ind, vals = self.hdg_wrapper.sparse_stiff_mat(0.)
    self.stiff_full = sp.sparse.csr_matrix( (vals, (row_ind, col_ind)),
                                          shape=(self.long_vector_size(), self.long_vector_size()))
    col_ind = np.array([self.index_dict.get(ind, -1) for ind in col_ind])
    row_ind = np.array([self.index_dict.get(ind, -1) for ind in row_ind])
    valid_mask = (col_ind > -1) & (row_ind > -1)
    col_ind = col_ind[valid_mask]
    row_ind = row_ind[valid_mask]
    vals = -np.array(vals)[valid_mask]
    return sp.sparse.csr_matrix((vals, (row_ind, col_ind)), shape=(system_size, system_size))
  def mass_mat(self, sigma = 1.):
    system_size = self.short_vector_size()
    _, _, vals_sig = self.hdg_wrapper.sparse_stiff_mat(sigma)
    col_ind, row_ind, vals = self.hdg_wrapper.sparse_stiff_mat(0.)

    # Vectorized index lookup with np.vectorize
    col_ind = np.array([self.index_dict.get(ind, -1) for ind in col_ind])
    row_ind = np.array([self.index_dict.get(ind, -1) for ind in row_ind])

    # Create a mask for valid indices
    valid_mask = (col_ind > -1) & (row_ind > -1)

    # Apply mask to filter arrays
    col_ind = col_ind[valid_mask]
    row_ind = row_ind[valid_mask]
    vals = np.array(vals)[valid_mask]
    vals_sig = np.array(vals_sig)[valid_mask]

    # Vectorized computation
    vals = (vals - vals_sig) / -sigma

    # Create the sparse matrix
    return sp.sparse.csr_matrix((vals, (row_ind, col_ind)), shape=(system_size, system_size))


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx_SI(poly_degree, dimension, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")

  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "NonlinearEigenvalue"
  const.topology        = "File<1,3>"
  const.geometry        = "File<1,3>"
  const.node_descriptor = "File<1,3>"
  const.cython_replacements = ["string", "string"]
  const.local_solver    = "TimoshenkoBeamEigs<1,3" + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",LocalSolver::TimoschenkoBeamParametersDefault,double>"
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/fiber_network_1000.geo" )
  HDG_wrapper.refine(2 ** iteration);

  helper = helper_ev_approx(HDG_wrapper)
  
  system_size = helper.short_vector_size()
  # Stiff       = LinearOperator( (system_size,system_size), matvec= helper.multiply_stiff )
  Stiff       = helper.stiff_mat()
  # Mass        = LinearOperator( (system_size,system_size), matvec= helper.multiply_mass )
  Mass        = helper.mass_mat()
  # Mass_inv    = LinearOperator( (system_size,system_size), matvec= helper.inv_mass )
  # ShiftedInv  = LinearOperator( (system_size,system_size), matvec= helper.shifted_inverse )



  points = np.loadtxt("domains/fiber_network_1000_points.txt")
  helper_precond = HyperHDG.fiber_network.precond(points, [2**3, 2**3], repeat=6)
  def precond_mult( vec_x ):
    result = helper_precond.precond(helper.stiff_full, helper.long_vector(vec_x))
    return helper.short_vector(result)
  B = sp.sparse.linalg.LinearOperator( (system_size,system_size), matvec= precond_mult )

  def mult_a_inv(vector):
    solution, n_iter = sp.sparse.linalg.cg(Stiff, vector, rtol=1e-10, M=B)
    return solution

  Stiff_inv = LinearOperator( (system_size,system_size), matvec= mult_a_inv )

  def mult_m_inv(vector):
    print("Eval m_inv")
    solution, n_iter = sp.sparse.linalg.cg(Mass, vector, rtol=1e-10)
    return solution

  Mass_inv = LinearOperator( (system_size,system_size), matvec= mult_m_inv )


  end_time = datetime.now()
  print("Setup done at", end_time, "after", end_time-start_time)

  # vals, vecs = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, sigma= sigma, which='LM', OPinv= ShiftedInv)
  # vals, vecs = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, which='SM')
  vals, vecs = sp_lin_alg.eigsh(Stiff_inv, k=1, M=Mass_inv, which='LM', Minv=Mass)
  vals[0] = 1/vals[0]

  end_time = datetime.now()
  print("Calculation done at", end_time, "after", end_time-start_time)

  error = np.absolute(vals[0])
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_hypergraph_convergence_eigenvalue_shifted_inverse.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  solution = helper.long_vector([x[0].real for x in vecs])
  
  HDG_wrapper.plot_option( "fileName" , "diff_eig_shifi-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(solution, vals[0].real)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  return vals[0].real, solution, error
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for iteration in range(10):
      eigenvalue_approx_SI(poly_degree, 1, iteration, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
