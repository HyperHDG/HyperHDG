from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration, refinement, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Elliptic"
  const.topology        = "Cubic<" + str(dimension) + ",3>"
  const.geometry        = "UnitCube<" + str(dimension) + ",3,double>"
  const.node_descriptor = "Cubic<" + str(dimension) + ",3>"
  const.local_solver    = "Diffusion<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(6) + ",HG<" + str(dimension) + ">::TestParametersQuadEllipt,double>"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( [2 ** iteration] * 3 )
  
  HDG_wrapper.refine( 2 ** refinement )

  vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1.)

  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
  if num_iter != 0:
    raise RuntimeError("Linear solver did not converge!")

  error = HDG_wrapper.errors(vectorSolution)[0]
  print("Refinement: ", refinement, " Error: ", error)
  f = open("output/diffusion_convergence_hypergraph_elliptic.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Refinement = " + str(refinement) + ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "diff_hyg_conv-" + str(dimension) + "-" + str(refinement) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  iteration = 2
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,4):
      print("Dimension is ", dimension, "\n")
      for refinement in range(4):
        try:
          diffusion_test(poly_degree, dimension, iteration, refinement, debug_mode)
        except RuntimeError as error:
          print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
