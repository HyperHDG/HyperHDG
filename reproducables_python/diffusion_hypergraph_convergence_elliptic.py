# Python running example for HDG solution of an elasticity problem on a superaggregate!

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
sys.path.append(os.path.dirname(__file__) + "/..")


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration, debug_mode=False):
  # Print starting time of diffusion test.
  start_time = datetime.now()
  print("Starting time is", start_time)

  # Predefine problem to be solved.
  problem = "GlobalLoop::Elliptic < Topology::Cubic<" + str(dimension) + ",3>, " \
          + "Geometry::UnitCube<" + str(dimension) + ",3,double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + ",3>, " \
          + "LocalSolver::Diffusion<" + str(dimension) + "," + str(poly_degree) + "," \
          + str(2*poly_degree) + ",HG<" + str(dimension) + ">::TestParametersQuadEllipt,double> >"
  filenames = [ "HyperHDG/geometry/unit_cube.hxx" , "HyperHDG/node_descriptor/cubic.hxx", \
                "HyperHDG/local_solver/diffusion_ldgh.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["elliptic_loop", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames, \
           debug_mode )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * 3 )

  # Generate right-hand side vector.
  vectorRHS = np.multiply(HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector()), -1.)

  # Define LinearOperator in terms of C++ functions to use scipy in a matrix-free fashion.
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
  if num_iter != 0:
    raise RuntimeError("Linear solver did not converge!")

  # Print error.
  error = HDG_wrapper.calculate_L2_error(vectorSolution)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_hypergraph_convergence_elliptic.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "diff_conv_hyg-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  
  # Print ending time of diffusion test.
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,4):
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
