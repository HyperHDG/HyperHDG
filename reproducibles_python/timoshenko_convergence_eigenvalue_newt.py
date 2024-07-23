from __future__ import print_function

import numpy as np
import scipy as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator
import random

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_newt():
  def __init__(self, hdg_wrapper):
    self.hdg_wrapper = hdg_wrapper
    dirichlet_indices = hdg_wrapper.dirichlet_nodes()
    self.dirichlet_indices = [index * 6 + comp for index in dirichlet_indices for comp in range(6)]
    self.index_vector      = [-1] * self.short_vector_size()
    self.vector_index      = [-1] * self.long_vector_size()
    self.val = [0.] * self.short_vector_size()
    n_indices = 0
    for i in range(len(self.index_vector)):
      while n_indices < len(self.dirichlet_indices) and \
            i + n_indices == self.dirichlet_indices[n_indices]:
        n_indices = n_indices + 1
      self.index_vector[i] = i + n_indices
      self.vector_index[i + n_indices] = i
    assert -1 not in self.index_vector
    self.index_dict = {value: idx for idx, value in enumerate(self.index_vector)}
  def set_val(self, val):
    self.val         = val
  def long_vector_size(self):
    return self.hdg_wrapper.size_of_system() + 1
  def short_vector_size(self):
    return self.hdg_wrapper.size_of_system() - len(self.dirichlet_indices) + 1
  def long_vector(self, vector):
    return [vector[x] if x > -1 else 0. for x in self.vector_index]
  def short_vector(self, vector):
    return [vector[x] for x in self.index_vector]
  def eval_residual(self, vector):
    vec = self.short_vector(self.hdg_wrapper.trace_to_flux(self.long_vector(vector), vector[-1]))
    # vec[-1] = np.linalg.norm(vector[:-1]) ** 2 - 1.
    vec[-1] = 0.
    return vec
  def eval_jacobian(self, vector):
    vec     = self.short_vector(self.hdg_wrapper.jacobian_of_trace_to_flux \
             ( self.long_vector(vector), vector[-1], self.long_vector(self.val), self.val[-1] ))
    vec[-1] = 2. * np.dot(vector[:-1],self.val[:-1])
    return vec
  def jacobian_mat(self, vector):
    system_size = self.short_vector_size()
    col_ind, row_ind, vals = self.hdg_wrapper.sparse_jacobi_mat(self.long_vector(vector))
    col_ind = np.array([self.index_dict.get(ind, -1) for ind in col_ind])
    row_ind = np.array([self.index_dict.get(ind, -1) for ind in row_ind])
    valid_mask = (col_ind > -1) & (row_ind > -1)
    col_ind = col_ind[valid_mask]
    row_ind = row_ind[valid_mask]
    vals = np.array(vals)[valid_mask]
    valid_mask = (col_ind == self.short_vector_size()-1)
    col_ind = np.concatenate((col_ind, row_ind[valid_mask]))
    vals = np.concatenate((vals, vals[valid_mask]))
    row_ind = np.concatenate((row_ind, np.array(([self.short_vector_size()-1] * np.sum(valid_mask)))))
    matrix = sp.sparse.csr_matrix((vals, (row_ind, col_ind)), shape=(system_size, system_size))
    return matrix



  # def eval_residual(self, vector):
  #   vec  = self.hdg_wrapper.trace_to_flux(vector, vector[len(vector)-1])
  #   temp = vector[len(vector)-1]
  #   vector[len(vector)-1] = 0.
  #   vec[len(vec)-1]       = np.linalg.norm(vector) ** 2 - 1.
  #   vector[len(vector)-1] = temp
  #   return vec
  # def eval_jacobian(self, vector):
  #   vec  = self.hdg_wrapper.jacobian_of_trace_to_flux \
  #            ( vector, vector[len(vector)-1], self.val, self.val[len(self.val)-1] )
  #   temp = self.val[len(self.val)-1]
  #   # self.val[-1] = vector[-1]
  #   # vec = np.add(vec, self.eval_residual(self.val))
  #   self.val[len(self.val)-1] = 0.
  #   vec[len(vec)-1]           = 2. * np.dot(vector,self.val)
  #   self.val[len(self.val)-1] = temp
  #   return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_newt(poly_degree, dimension, iteration, initial="default", debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")

  n_newton_steps = 10 ** 2
  scaling_fac    = 1e-2
  norm_exact     = 1e-8
  alpha          = 0.1
  beta           = 0.5

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
  # const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/cross.geo" )
  HDG_wrapper.refine(2 ** iteration);
  
  helper = helper_ev_newt(HDG_wrapper)

  system_size = helper.short_vector_size()
  A = LinearOperator( (system_size,system_size), matvec= helper.eval_jacobian )


  
  if initial == "default":
    vectorSolution = helper.long_vector([0.] * system_size)
    vectorSolution = HDG_wrapper.make_initial(vectorSolution)
    vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    vectorSolution[-1] = dimension * (np.pi ** 2) + 0. * random.randint(-100,100)
  else:
    vectorSolution = initial
    # temp = initial[len(vectorSolution)-1]
    # vectorSolution[len(vectorSolution)-1] = 0.
    # vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    # vectorSolution[len(vectorSolution)-1] = temp
    # vectorSolution = np.array([entry + 1e-6 * random.randint(-100,100) for entry in vectorSolution])

  # print(vectorSolution)

  vectorSolution = helper.short_vector(vectorSolution)
  # print(vectorSolution)
  # A = helper.jacobian_mat(vectorSolution)
  residual = helper.eval_residual( vectorSolution )
  norm_res = np.linalg.norm( residual )

  # print(vectorSolution)

  # solution = sp.optimize.newton(helper.eval_residual, vectorSolution, rtol=1e-10, full_output=True)
  # vectorSolution = solution[0]



  for newton_step in range(n_newton_steps):
    
    # helper.set_val(vectorSolution)
    norm_old = norm_res
    gamma = 1

    print(norm_old)

    A = helper.jacobian_mat(vectorSolution).todense()

    [vectorUpdate, num_iter] = sp_lin_alg.cg(A, residual, rtol= scaling_fac * norm_res)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations in step", newton_step ," with norm bound", scaling_fac * norm_res, ".")
      [vectorUpdate, num_iter] = sp_lin_alg.bicgstab(A, residual, rtol= scaling_fac * norm_res)
      if num_iter != 0:
        print("BiCGStab also failed with ", num_iter, "iterations")
        raise RuntimeError("Linear solvers did not converge!")
    
    vectorHelper = np.subtract(vectorSolution, vectorUpdate)
    residual = helper.eval_residual(vectorHelper)
    norm_res = np.linalg.norm( residual )
    
    while norm_res > (1 - alpha * gamma) * norm_old:
      if gamma < 1e-6:
        raise RuntimeError("Newton step is too small!")
      gamma = beta * gamma
      vectorHelper = np.subtract(vectorSolution, np.multiply(vectorUpdate, gamma))
      residual = helper.eval_residual(vectorHelper)
      norm_res = np.linalg.norm( residual )
    
    vectorSolution = vectorHelper
    vectorSolution[:-1] = np.divide(vectorSolution[:-1], np.linalg.norm(vectorSolution[:-1]))
    residual = helper.eval_residual(vectorHelper)
    norm_res = np.linalg.norm( residual )
    
    if norm_res < norm_exact:
      print("Newton solver converged after ", newton_step+1, " iterations with residual norm", norm_res, ".")
      break
      
  if norm_res >= norm_exact:
    raise RuntimeError("Newton solver did not converge!")
  
  error = np.absolute(vectorSolution[-1] - np.pi**2 / 4)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_convergence_eigenvalue_newton.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "diff_eig_newt-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(helper.long_vector(vectorSolution),vectorSolution[-1])
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  return vectorSolution[-1], vectorSolution, error
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,2):
    print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
    for iteration in range(1):
      try:
        eigenvalue_newt(poly_degree, 1, iteration, "default", debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
