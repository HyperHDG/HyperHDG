# Python running example for HDG solution of an parabolic diffusion equation with implicit Euler!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Import package to print date and time of start and end of program.
from datetime import datetime

# Correct the python paths!
import os, sys

# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(theta, poly_degree, refinement, debug_mode=False):
  # Print starting time of diffusion test.
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  # Config time stepping.
  time_steps      = 10 ** 3
  time_end        = 5
  delta_time      = time_end / time_steps
  output_interval = 20
  
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
    ",LeVequeHG,double>"
  const.cython_replacements = ["string", "string", "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/advection.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/leVeque_hg.geo", lsol_constr= [0.,theta,delta_time] )
  HDG_wrapper.refine( 2 ** refinement )

  # Generate right-hand side vector.
  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  
  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "leVeque_hyg" + str(theta) + "-" + str(poly_degree) \
    + "-" + str(refinement) )
  HDG_wrapper.plot_option( "printFileNumber" , "true" )
  # HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time_end)

  # For loop over the respective time-steps.
  for time_step in range(time_steps):
    
    # Assemble right-hand side vextor and "mass_matrix * old solution".
    vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
                 (time_step+1) * delta_time), -1.)

    # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
    # [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
    # if num_iter != 0:
    #   print("CG failed with a total number of ", num_iter, " iterations in time step ", time_step, \
    #         ". Trying GMRES!")
    [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-13)
    if num_iter != 0:
      # print("GMRES also failed with a total number of ", num_iter, "iterations.")
      [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,tol=1e-13)
      if num_iter != 0:
        print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
        raise RuntimeError("All linear solvers did not converge!")

    HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)

    if (time_step+1) % output_interval == 0:
      HDG_wrapper.plot_solution(vectorSolution, time_end)
    
  # Print error.
  error = HDG_wrapper.errors(vectorSolution, time_end)[0]
  print( "Iteration: ", refinement, " Error: ", error )
  f = open("output/advection_convergence_rotation_theta"+str(theta)+".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Theta = " + str(theta) \
          + ". Iteration = " + str(refinement) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  # HDG_wrapper.plot_option( "fileName" , "leVeque_hyg" + str(theta) + "-" + str(poly_degree) \
  #   + "-" + str(refinement) )
  # HDG_wrapper.plot_option( "printFileNumber" , "true" )
  # HDG_wrapper.plot_option( "scale" , "0.95" )
  # HDG_wrapper.plot_solution(vectorSolution, time_end)
  
  # Print ending time of diffusion test.
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for theta in [0.5, 1.]:
    print("\n Theta is set to be ", theta, "\n\n")
    for poly_degree in range(3):
      print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
      for refinement in range(5,8):
        try:
          diffusion_test(theta, poly_degree, refinement, debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
