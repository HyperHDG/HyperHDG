from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys

## NOTE:
# Running the ..._approx.py -files currently produce an error:
# - vals, vecs = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, sigma= sigma, which='LM', OPinv = ShiftedInv)
# - - arpack.py, line 1700, in eigsh params.iterate()
# - - arpack.py, line 573, in iterate raise ArpackError(self.info, infodict=self.iterate_infodict)
# - - - scipy.sparse.linalg._eigen.arpack.arpack.ArpackError: ARPACK error -9: Starting vector is zero.

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
    vec = np.multiply(self.hdg_wrapper.trace_to_flux( self.long_vector(vector) ), -1.)
    vec = self.short_vector(vec)
    return vec
  def multiply_mass(self, vector):
    vec = self.short_vector( self.hdg_wrapper.trace_to_mass_flux( self.long_vector(vector) ) )
    return vec
  def shifted_mult(self, vector):
    vec = np.multiply( self.multiply_mass(vector), self.sigma)
    vec = np.subtract( self.multiply_stiff(vector), vec)
    return vec
  def shifted_inverse(self, vector):
    if not hasattr(self, 'ShiftOp'):
      system_size  = self.short_vector_size()
      self.ShiftOp = LinearOperator( (system_size,system_size), matvec=self.shifted_mult )
    [vec, n_it] = sp_lin_alg.cg(self.ShiftOp, vector, atol=1e-9)
    assert n_it == 0
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx_MA(poly_degree, dimension, iteration, level, debug_mode=False):
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

  try:
    from star_shape_functions import get_star_domain
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../domains")
    from star_shape_functions import get_star_domain

  # Determine the width of the star, and if level = 0 -> limiting case
  if (level == 0):
    width = 1
    hy_dim = 1
    filename = f"star_{dimension}_lim.geo"
  else:
    width = 2 ** level + 1
    hy_dim = dimension
    filename = get_star_domain(dimension, level)

  const                 = HyperHDG.config()
  const.global_loop     = "MassApproxEigenvalue"
  const.topology        = "File<" + str(hy_dim) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.geometry        = "File<" + str(hy_dim) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.node_descriptor = "File<" + str(hy_dim) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.local_solver    = "DiffusionEigs<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",ThicknessEigs<" + str(width) + ">::TestParametersEigs,double>"
  const.cython_replacements = ["string", "string"]
  const.include_files   = ["reproducibles_python/parameters/diffusion_star.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP(os.path.dirname(os.path.abspath(__file__)) + "/../domains/" + filename)
  HDG_wrapper.refine(2 ** iteration)
  
  helper = helper_ev_approx(HDG_wrapper, sigma)
  
  system_size = helper.short_vector_size()

  Stiff       = LinearOperator( (system_size,system_size), matvec= helper.multiply_stiff )
  Mass        = LinearOperator( (system_size,system_size), matvec= helper.multiply_mass )
  ShiftedInv  = LinearOperator( (system_size,system_size), matvec= helper.shifted_inverse )

  vals, vecs = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, sigma= sigma, which='LM', OPinv = ShiftedInv)

  eig_val = vals[0]
  print("Iteration: ", iteration, " Eigenvalue: ", eig_val)
  f = open("output/diffusion_star_eigenvalue_approx.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". Eigenvalue = " + str(eig_val) + ".\n")
  f.close()
  
  solution = helper.long_vector([x[0].real for x in vecs])
  
  HDG_wrapper.plot_option( "fileName" , "diff_star_apprx_dim-" + str(dimension) + "-lvl-" + str(level) + "-" + str(iteration))
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(solution);
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  return vals[0].real, solution, eig_val
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1, 3):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(2, 3):
      print("\n Dimension is set to be ", dimension, "\n\n")
      for level in range(1, 4):
        print("Level is ", level, "\n")
        for iteration in range(2, 5):
          eigenvalue_approx_MA(poly_degree, dimension, iteration, level, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
