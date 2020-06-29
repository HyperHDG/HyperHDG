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


class helper_class():
  def __init__(self, hdg_class, delta_time):
    self.hdg_wrapper = hdg_class
    self.delta_time = delta_time
  def multiply(self, vec):
    vec1 = np.multiply(self.hdg_wrapper.matrix_vector_multiply(vec), -self.delta_time)
    vec2 = self.hdg_wrapper.mass_matrix_multiply(vec)
    return np.add(vec1, vec2)

# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(dimension, iteration):
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Diffusion<" + str(dimension) + ",1,2,TestParametersSinParab,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/Diffusion.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames, True )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )
  time_steps  = 40 * iteration
  delta_time  = 1 / time_steps

  helper = helper_class(HDG_wrapper, delta_time)

  # Generate right-hand side vector.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.mass_matrix_multiply )
  vectorSolution = HDG_wrapper.initial_flux_vector(HDG_wrapper.return_zero_vector())
  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorSolution, tol=1e-9)
  
  # Generate right-hand side vector.
  vectorRHS = [-i for i in HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector())]
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )
  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-9)

  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= helper.multiply )

  print("Delta time:" + str(delta_time))

  for time_step in range(time_steps):

    vectorRHS = \
      np.multiply(HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector()), delta_time)
    vectorSolution = HDG_wrapper.mass_matrix_multiply(vectorSolution)

    # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
    [vectorSolution, num_iter] = sp_lin_alg.cg(A, np.add(vectorRHS,vectorSolution), tol=1e-9)
    
    if num_iter != 0:
      print("The linear solver (conjugate gradients) failed with a total number of ",
            num_iter, " iterations.")

  # Print error.
  print("Error: " + str(HDG_wrapper.calculate_L2_error(vectorSolution, 0.)))
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "diff_c-" + str(dimension) + "-" + str(iteration) );
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(vectorSolution);
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for dimension in range(1,4):
    for iteration in range(1, 8 - dimension):
      diffusion_test(dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()