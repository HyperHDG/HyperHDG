from __future__ import print_function

import numpy as np
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
    self.val         = [0] * (hdg_wrapper.size_of_system() + 1)
  def set_val(self, val):
    self.val         = val
  def eval_residual(self, vector):
    vec  = self.hdg_wrapper.trace_to_flux(vector, vector[len(vector)-1])
    temp = vector[len(vector)-1]
    vector[len(vector)-1] = 0.
    vec[len(vec)-1]       = np.linalg.norm(vector) ** 2 - 1.
    vector[len(vector)-1] = temp
    return vec
  def eval_jacobian(self, vector):
    vec  = self.hdg_wrapper.jacobian_of_trace_to_flux \
             ( vector, vector[len(vector)-1], self.val, self.val[len(self.val)-1] )
    temp = self.val[len(self.val)-1]
    self.val[len(self.val)-1] = 0.
    vec[len(vec)-1]           = 2. * np.dot(vector,self.val)
    self.val[len(self.val)-1] = temp
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_newt(poly_degree, dimension, level, iteration, initial="default", debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")

  n_newton_steps = 10 ** 1
  scaling_fac    = 1e-2
  norm_exact     = 1e-9
  alpha          = 0.1
  beta           = 0.5

  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG

  try:
    from L_shape_functions import get_l_domain
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../domains")
    from L_shape_functions import get_l_domain

  if (level == 0):
    hy_dim = 1
    filename = "L_lim.geo"
  else:
    hy_dim = 2
    filename = get_l_domain(level)

  const                 = HyperHDG.config()
  const.global_loop     = "NonlinearEigenvalue"
  const.topology        = "File<" + str(hy_dim) + ",2,std::vector,Point<2,double> >"
  const.geometry        = "File<" + str(hy_dim) + ",2,std::vector,Point<2,double> >"
  const.node_descriptor = "File<" + str(hy_dim) + ",2,std::vector,Point<2,double> >"
  const.local_solver    = "DiffusionEigs<" + str(hy_dim) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",ThicknessEigs<"+str(level)+">::TestParametersEigs,double>"
  const.cython_replacements = ["string", "string"]
  const.include_files   = ["reproducibles_python/parameters/diffusion_L.hxx"]
  const.debug_mode      = debug_mode

  # filename = get_l_domain(level) # "L_domain_dir.geo"

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP(os.path.dirname(os.path.abspath(__file__)) + "/../domains/" + filename)
  HDG_wrapper.refine(2 ** iteration)
  
  helper = helper_ev_newt(HDG_wrapper)

  system_size = HDG_wrapper.size_of_system() + 1
  A = LinearOperator( (system_size,system_size), matvec= helper.eval_jacobian )
  
  if initial == "default":
    vectorSolution = [0] * system_size
    vectorSolution = HDG_wrapper.make_initial(vectorSolution)
    vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    vectorSolution[system_size-1] = dimension * (np.pi ** 2) + 1e-3 * random.randint(-100,100)
  else:
    vectorSolution = initial
    temp = initial[len(vectorSolution)-1]
    vectorSolution[len(vectorSolution)-1] = 0.
    vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    vectorSolution[len(vectorSolution)-1] = temp
  
  residual = helper.eval_residual( vectorSolution )
  norm_res = np.linalg.norm( residual )
  
  for newton_step in range(n_newton_steps):
    
    helper.set_val(vectorSolution)
    norm_old = norm_res
    gamma = 1
    
    [vectorUpdate, num_iter] = sp_lin_alg.gmres(A, residual, atol= scaling_fac * norm_res)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations.")
      [vectorUpdate, num_iter] = sp_lin_alg.bicgstab(A, residual, atol= scaling_fac * norm_res)
      if num_iter != 0:
        print("BiCGStab also failed with ", num_iter, "iterations")
        raise RuntimeError("Linear solvers did not converge!")
    
    vectorHelper = np.subtract(vectorSolution, vectorUpdate)
    residual = helper.eval_residual(vectorHelper)
    norm_res = np.linalg.norm( residual )
    
    while norm_res > (1 - alpha * gamma) * norm_old:
      if gamma < 1e-4:
        raise RuntimeError("Newton step is too small!")
      gamma = beta * gamma
      vectorHelper = np.subtract(vectorSolution, np.multiply(vectorUpdate, gamma))
      residual = helper.eval_residual(vectorHelper)
      norm_res = np.linalg.norm( residual )
    
    vectorSolution = vectorHelper
    
    if norm_res < norm_exact:
      print("Newton solver converged after ", newton_step+1, " iterations.")
      break
      
  if norm_res >= norm_exact:
    raise RuntimeError("Newton solver did not converge!")

  eig_val = vectorSolution[system_size-1]
  print("Iteration: ", iteration, " Eigenvalue: ", eig_val)
  f = open("output/diffusion_L_convergence_eigenvalue_newton.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". Eigenvalue = " + str(eig_val) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "diff_L_eig_newt_pd-" + str(poly_degree) + "-lvl-" + str(level) + "-" + str(iteration))
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution,vectorSolution[system_size-1])
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  return vectorSolution[system_size-1], vectorSolution, eig_val
  

# --------------------------------------------------------------------------------------------------
# Function main.diffusion_convergence_eigenvalue_newt
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1, 3):
    print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
    for level in range(0, 3):
      print("\nLevel is ", level, "\n")
      for iteration in range(2,4):
        try:
          eigenvalue_newt(poly_degree, 2, level, iteration, "default", debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
