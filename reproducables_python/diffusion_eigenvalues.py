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
def diffusion_test(dimension, iteration):
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Diffusion<" + str(dimension) + ",1,2,TestParametersHomo,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/Diffusion.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames, True )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )

  # Generate right-hand side vector.
  vectorRHS = [-i for i in HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector())]

  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )
  Mass = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.mass_matrix_multiply )

  print( A * [1, 0, 0] )
  print( A * [0, 1, 0] )
  print( A * [0, 0, 1] )
  print("Make mass")
  print( Mass * [1, 0, 0] )
  print( Mass * [0, 1, 0] )
  print( Mass * [0, 0, 1] )


  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
 # [vals, vecs] = sp_lin_alg.eigs(A, k=1, M=Mass, which='LR', tol=1e-2)

  # Print error.
 # print(vals)
 # print(vecs)
  
  # Plot obtained solution.
 # HDG_wrapper.plot_option( "fileName" , "diff_e-" + str(dimension) + "-" + str(iteration) );
 # HDG_wrapper.plot_option( "printFileNumber" , "false" );
 # HDG_wrapper.plot_option( "scale" , "0.95" );
 # HDG_wrapper.plot_solution([i.real for i in vecs]);
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for dimension in range(1,2):
    for iteration in range(1, 2):#6 - dimension):
      diffusion_test(dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
