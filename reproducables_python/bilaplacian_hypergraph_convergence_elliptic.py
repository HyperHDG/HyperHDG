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
def diffusion_test(poly_degree, dimension, iteration):
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + ",3>, " \
          + "Geometry::UnitCube<" + str(dimension) + ",3,double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + ",3>, " \
          + "bilaplacian<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersQuadEllipt" + str(dimension) + ",double> >"
  filenames = [ "HyperHDG/geometry/cubic.hxx" , "HyperHDG/node_descriptor/cubic.hxx", \
                "HyperHDG/local_solver/bilaplacian.hxx", \
                "reproducables_python/parameters/bilaplacian.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * 3 )

  # Generate right-hand side vector.
  vectorRHS = np.multiply(HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector()), -1.)

  # Define LinearOperator in terms of C++ functions to use scipy in a matrix-free fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)

  # Print error.
  error = HDG_wrapper.calculate_L2_error(vectorSolution)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/bilaplacian_hypergraph_convergence_elliptic.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "bil_conv_hyg-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,4):
      print("Dimension is ", dimension, "\n")
      for iteration in range(6):
        diffusion_test(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
