from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys

# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, level, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG

  try:
    from star_shape_functions import get_star_domain
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../domains")
    from star_shape_functions import get_star_domain

  # Determine the width of the star
  if (level == 0):
    width = 1
  else:
    width = 2**level + 1

  const                 = HyperHDG.config()
  const.global_loop     = "Elliptic"
  const.topology        = "File<" + str(dimension) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.geometry        = "File<" + str(dimension) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.node_descriptor = "File<" + str(dimension) + "," + str(dimension) + ",std::vector,Point<" + str(dimension) + ",double> >"
  const.local_solver    = "Diffusion<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",Thickness<" + str(width) + ">::TestParametersSinEllipt, double>"
  const.cython_replacements = ["string", "string"]
  const.include_files   = ["reproducibles_python/parameters/diffusion_star.hxx"]
  const.debug_mode      = debug_mode

  filename = get_star_domain(dimension, level)

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/" + filename)
  HDG_wrapper.refine( 2 ** iteration )

  vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
  if num_iter != 0:
    print("CG solver failed with a total number of ", num_iter, "iterations.")
    [vectorSolution, num_iter] = sp_lin_alg.gmres(A, vectorRHS, tol=1e-13)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations.")
      raise RuntimeError("Linear solvers did not converge!")

  error = HDG_wrapper.errors(vectorSolution)[0]
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_plus_level_" + str(level) + ".txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Iteration = " + str(iteration) + \
          ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "diff_star_dim_" + str(dimension) + "_level-" + str(level) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1, 4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1, 4):
      print("\n Dimension is ", dimension, "\n")
      for level in range(0, 4):
        print("\nLevel is ", level, "\n")
        for iteration in range(3):
          try:
            diffusion_test(poly_degree, dimension, level, iteration, debug_mode)
          except RuntimeError as error:
            print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
