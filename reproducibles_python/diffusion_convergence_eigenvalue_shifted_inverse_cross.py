from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_approx():
  def __init__(self, hdg_wrapper, sigma):
    self.hdg_wrapper       = hdg_wrapper
    self.sigma             = sigma
    self.dirichlet_inidces = hdg_wrapper.dirichlet_indices()
    self.index_vector      = [-1] * self.short_vector_size()
    self.vector_index      = [-1] * self.long_vector_size()
    n_indices = 0
    for i in range(len(self.index_vector)):
      while n_indices < len(self.dirichlet_inidces) and \
            i + n_indices == self.dirichlet_inidces[n_indices]:
        n_indices = n_indices + 1
      self.index_vector[i] = i + n_indices
      self.vector_index[i + n_indices] = i
    assert -1 not in self.index_vector
  def long_vector_size(self):
    return self.hdg_wrapper.size_of_system()
  def short_vector_size(self):
    return self.hdg_wrapper.size_of_system() - len(self.dirichlet_inidces)
  def long_vector(self, vector):
    return [vector[x] if x > -1 else 0. for x in self.vector_index]
  def short_vector(self, vector):
    return [vector[x] for x in self.index_vector]
  def multiply_stiff(self, vector):
    vec = self.hdg_wrapper.trace_to_flux(self.long_vector(vector))
    vec = np.multiply(self.short_vector(vec), -1.)
    return vec
  def shifted_mult(self, vector):
    vec = self.hdg_wrapper.trace_to_flux(self.long_vector(vector), self.sigma)
    vec = np.multiply(self.short_vector(vec), -1.)
    return vec
  def multiply_mass(self, vector):
    vec = self.hdg_wrapper.trace_to_flux(self.long_vector(vector), self.sigma)
    vec = np.subtract( self.hdg_wrapper.trace_to_flux(self.long_vector(vector)) , vec )
    vec = np.multiply( self.short_vector(vec) , -1. / self.sigma )
    return vec
  def shifted_inverse(self, vector):
    if not hasattr(self, 'ShiftOp'):
      system_size  = self.short_vector_size()
      self.ShiftOp = LinearOperator( (system_size,system_size), matvec=self.shifted_mult )
    [vec, n_it] = sp_lin_alg.cg(self.ShiftOp, vector, tol=1e-9)
    assert n_it == 0
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx_SI(poly_degree, dimension, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  exact_eigenval = dimension * (np.pi ** 2)
  sigma          = exact_eigenval / 2.

  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "ShiftedInverseEigenvalue"
  const.topology        = "File<1,2,std::vector,Point<2,double> >"
  const.geometry        = "File<1,2,std::vector,Point<2,double> >"
  const.node_descriptor = "File<1,2,std::vector,Point<2,double> >"
  const.local_solver    = "DiffusionEigs<1," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",TestParametersEigs,double>"
  const.cython_replacements = ["string", "string"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/cross.geo" )
  HDG_wrapper.refine( 2 ** iteration )

  
  helper = helper_ev_approx(HDG_wrapper, sigma)
  
  system_size = helper.short_vector_size()
  Stiff       = LinearOperator( (system_size,system_size), matvec= helper.multiply_stiff )
  Mass        = LinearOperator( (system_size,system_size), matvec= helper.multiply_mass )
  ShiftedInv  = LinearOperator( (system_size,system_size), matvec= helper.shifted_inverse )
    
  vals, vecs = sp_lin_alg.eigsh(Stiff, k=2, M=Mass, sigma= sigma, which='LM', OPinv= ShiftedInv)

  error = vals[0]
  print("Iteration: ", iteration, " EV1: ", error, " EV2: ", vals[1])
  f = open("output/diffusion_convergence_eigenvalue_shifted_inverse.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". EV2 = " + str(error) + ". EV2 = " + str(vals[1]) + "\n")
  f.close()
  
  solution = helper.long_vector([x[1].real for x in vecs])
  
  HDG_wrapper.plot_option( "fileName" , "diff_eig_shifi-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(solution, vals[1].real)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  return vals[0].real, solution, error
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(2,6):
        eigenvalue_approx_SI(poly_degree, dimension, iteration, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
