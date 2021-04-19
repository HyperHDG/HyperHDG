from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  theta       = 1.
  time_steps  = 10 ** 4
  delta_time  = 1 / time_steps
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Parabolic"
  const.topology        = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
  const.geometry        = "UnitCube<" + str(dimension) + "," + str(dimension) + ",double>"
  const.node_descriptor = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
  const.local_solver    = "DiffusionParab<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",TestParametersSinParab,double>"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( [2 ** iteration] * dimension, lsol_constr= [1.,theta,delta_time] )

  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  for time_step in range(time_steps):
    
    vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
                 (time_step+1) * delta_time), -1.)
    
    [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
    if num_iter != 0:
      print("CG failed with a total number of ", num_iter, " iterations in time step ", time_step, \
            ". Trying GMRES!")
      [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-13)
      if num_iter != 0:
        print("GMRES also failed with a total number of ", num_iter, "iterations.")
        [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,tol=1e-13)
        if num_iter != 0:
          print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
          raise RuntimeError("All linear solvers did not converge!")

    HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)
    
  error = HDG_wrapper.errors(vectorSolution, 1.)[0]
  print( "Iteration: ", iteration, " Error: ", error )
  f = open("output/diffusion_convergence_parabolic_theta"+str(theta)+".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "diff_conv_parab" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, 1.)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(6):
        try:
          diffusion_test(poly_degree, dimension, iteration, debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
