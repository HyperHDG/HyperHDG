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
def diffusion_test(poly_degree, dimension, iteration, debug_mode=False):
  # Print starting time of diffusion test.
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  # Config time stepping.
  theta       = 1.
  time_steps  = 10 ** 4
  time_end    = 1.
  delta_time  = time_end / time_steps
  
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
  const.local_solver    = "AdvectionParab<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",TestParametersSinParab,double>"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", \
    "double", "vector[double]"]
  const.include_files   = ["reproducables_python/parameters/advection.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP([2 ** iteration] * dimension, lsol_constr= [1.,theta,delta_time] )

  # Generate right-hand side vector.
  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  
  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

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
    
  # Print error.
  error = HDG_wrapper.errors(vectorSolution, time_end)[0]
  print( "Iteration: ", iteration, " Error: ", error )
  f = open("output/advection_convergence_parabolic_theta"+str(theta)+".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "adv_conv_parab" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time_end)
  
  # Print ending time of diffusion test.
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,2):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,2):
      print("Dimension is ", dimension, "\n")
      for iteration in range(0,8):
        try:
          diffusion_test(poly_degree, dimension, iteration, debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
