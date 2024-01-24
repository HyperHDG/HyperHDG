from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys

# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(theta, poly_degree, refinement, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  time_steps      = 10 ** 3
  time_end        = 5
  delta_time      = time_end / time_steps
  output_interval = 1
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Parabolic"
  const.topology        = "File<2,3,std::vector,Point<3,double> >"
  const.geometry        = "File<2,3,std::vector,Point<3,double> >"
  const.node_descriptor = "File<2,3,std::vector,Point<3,double> >"
  const.local_solver    = "AdvectionParab<2," + str(poly_degree) + "," + str(2*poly_degree) + \
    ",HyperGraphMovie,double>"
  const.cython_replacements = ["string", "string", "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/advection.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/leVeque_hg.geo", lsol_constr= [0.,theta,delta_time] )
  HDG_wrapper.refine( 2 ** refinement )

  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  HDG_wrapper.plot_option( "fileName" , "leVeque_hyg" + str(theta) + "-" + str(poly_degree) \
    + "-" + str(refinement) )
  HDG_wrapper.plot_option( "printFileNumber" , "true" )
  # HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time_end)

  for time_step in range(time_steps):
    
    vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
                 (time_step+1) * delta_time), -1.)

    [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-8)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations.")
      [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,tol=1e-8)
      if num_iter != 0:
        print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
        raise RuntimeError("All linear solvers did not converge!")

    HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)

    if (time_step+1) % output_interval == 0:
      HDG_wrapper.plot_solution(vectorSolution, time_end)
    
  error = HDG_wrapper.errors(vectorSolution, time_end)[0]
  print( "Iteration: ", refinement, " Error: ", error )
  f = open("output/advection_convergence_rotation_theta"+str(theta)+".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Theta = " + str(theta) \
          + ". Iteration = " + str(refinement) + ". Error = " + str(error) + ".\n")
  f.close()
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for theta in [0.5, 1.]:
    print("\n Theta is set to be ", theta, "\n\n")
    for poly_degree in range(1,4):
      print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
      for refinement in range(5,6):
        try:
          diffusion_test(theta, poly_degree, refinement, debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
