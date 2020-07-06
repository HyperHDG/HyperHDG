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
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_class():
  def __init__(self, hdg_wrapper, delta_time):
    self.hdg_wrapper = hdg_wrapper
    self.delta_time = delta_time
  def multiply(self, vec):
    vec_stiff = np.multiply(self.hdg_wrapper.matrix_vector_multiply(vec), -self.delta_time)
    vec_mass = self.hdg_wrapper.mass_matrix_multiply(vec)
    return np.add(vec_stiff, vec_mass)


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def bilaplacian_test(dimension, iteration):
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "bilaplacian<" + str(dimension) + ",1,2,TestParametersSinParab,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/bilaplacian.hxx", \
                "reproducables_python/parameters/bilaplacian.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames, True )
  
  # Config time stepping.
  time_steps  = 4 * (iteration+1) * (iteration+1)
  delta_time  = 1 / time_steps
  
  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension, tau= (2**iteration) ) # Why is this so good?
  helper = helper_class(HDG_wrapper, delta_time)

  # Generate right-hand side vector.
  vectorSolution = HDG_wrapper.initial_flux_vector(HDG_wrapper.return_zero_vector())
  
  
  # Generate right-hand side vector.
  vectorRHS = np.multiply(HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector()), -1.)
  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )
  # Solve "A * x = b" in matrix-free fashion using scipy's BiCGStab algorithm (much faster than CG).
  [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A, vectorRHS, tol=1e-14)
  
  
  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= helper.multiply )
  
  # For loop over the respective time-steps.
  for time_step in range(time_steps):
    
    # Assemble right-hand side vextor and "mass_matrix * old solution".
    vectorRHS = np.multiply( \
      HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector(), (time_step+1) * delta_time), \
      delta_time )
    vectorSolution = HDG_wrapper.mass_matrix_multiply(vectorSolution)

    # Solve "A * x = b" in matrix-free fashion using scipy's BiCGStab algorithm.
    [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A, np.add(vectorRHS,vectorSolution), tol=1e-14)
    if num_iter != 0:
      print("The linear solver (conjugate gradients) failed with a total number of ",
            num_iter, " iterations.")

  # Print error.
  print("Error: " + str(HDG_wrapper.calculate_L2_error(vectorSolution, 1.)))
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "bilap_c-" + str(dimension) + "-" + str(iteration) );
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(vectorSolution);
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for dimension in range(1,2):
    for iteration in range(1, 10 - dimension):
      bilaplacian_test(dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
