# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Correct the python paths!
import os, sys
sys.path.append(os.path.dirname(__file__) + "/..")


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def bilaplacian_test(poly_degree, dimension, iteration):
  
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "bilaplacian_parab<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersSinParab,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/bilaplacian.hxx", \
                "reproducables_python/parameters/bilaplacian.hxx" ]

  # Config time stepping.
  time_steps  = 10 ** 4
  delta_time  = 1 / time_steps

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]", \
            "double", "vector[double]"], filenames )
  
  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension, lsol_constr= [1.,1.,delta_time] )

  # Generate right-hand side vector.
  vectorSolution = HDG_wrapper.initial_flux_vector(HDG_wrapper.return_zero_vector())
  
  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )
  
  # For loop over the respective time-steps.
  for time_step in range(time_steps):
    
    if time_step > 0:
      HDG_wrapper.set_data(vectorSolution, time_step * delta_time)
    
    # Assemble right-hand side vextor and "mass_matrix * old solution".
    vectorRHS = np.multiply(HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector(), \
                 (time_step+1) * delta_time), -1.)
    
    # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
    [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
    if num_iter != 0:
      print("CG failed with a total number of ", num_iter, " iterations in time step ", time_step, \
            ". Trying GMRES!")
      [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-13)
      if num_iter != 0:
        print("GMRES also failed with a total number of ", num_iter, "iterations.")
        sys.exit("Program failed!")

  # Print error.
  error = HDG_wrapper.calculate_L2_error(vectorSolution, 1.)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/bilaplacian_convergence_parabolic.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "bilap_c_p-" + str(dimension) + "-" + str(iteration) );
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(vectorSolution, 1.);
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(6):
        bilaplacian_test(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
