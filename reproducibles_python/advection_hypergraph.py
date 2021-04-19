from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(tau, theta, poly_degree, edge_dim, space_dim, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")

  time_steps  = 10 ** 3
  time_end    = 10
  delta_time  = time_end / time_steps
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Parabolic"
  const.topology        = "File<" + str(edge_dim) + "," + str(space_dim) + ",std::vector,Point<" \
    + str(space_dim) + ",double> >"
  const.geometry        = "File<" + str(edge_dim) + "," + str(space_dim) + ",std::vector,Point<" \
    + str(space_dim) + ",double> >"
  const.node_descriptor = "File<" + str(edge_dim) + "," + str(space_dim) + ",std::vector,Point<" \
    + str(space_dim) + ",double> >"
  const.local_solver    = "AdvectionParab<" + str(edge_dim) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",injection,double>"
  const.cython_replacements = ["string", "string", "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/advection.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/injection_test.geo", lsol_constr= [tau,theta,delta_time] )

  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())

  HDG_wrapper.plot_option( "fileName" , "adv_injection" + str(tau) )
  HDG_wrapper.plot_option( "printFileNumber" , "true" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time_end)

  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  for time_step in range(time_steps):
    
    vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
                 (time_step+1) * delta_time), -1.)

    [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-13)
    if num_iter != 0:
      # print("GMRES also failed with a total number of ", num_iter, "iterations.")
      [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,tol=1e-13)
      if num_iter != 0:
        print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
        raise RuntimeError("All linear solvers did not converge!")

    HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)

    HDG_wrapper.plot_option( "fileName" , "adv_injection" + str(tau) )
    HDG_wrapper.plot_option( "printFileNumber" , "true" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, time_end)
    
  error = HDG_wrapper.errors(vectorSolution, time_end)[0]
  print( "Iteration: ", iteration, " Error: ", error )
  f = open("output/advection_injection_theta"+str(theta)+".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". edge_dim = " + str(edge_dim) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  edge_dim = 1
  space_dim = 2
  iteration = 2
  poly_degree = 3
  theta = 0.5
  for tau in [0, 1]:
    try:
      diffusion_test(tau, theta, poly_degree, edge_dim, space_dim, iteration, debug_mode)
    except RuntimeError as error:
      print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
