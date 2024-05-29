from __future__ import print_function

import numpy as np
import scipy.sparse as sp
# import scipy.sparse.linalg as sp_lin_alg
# from scipy.sparse.linalg import LinearOperator

# import matplotlib.pyplot as plt

from datetime import datetime

import os, sys

# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration, debug_mode=False):
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
  const.local_solver    = "TimoshenkoBeam<1,3,4,8>"
  const.topology        = "File<1,3>"
  const.geometry        = "File<1,3>"
  const.node_descriptor = "File<1,3>"
  const.cython_replacements = ["string", "string"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/cross.geo" )
  HDG_wrapper.refine(2 ** iteration);

  vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )
  
  system_size = HDG_wrapper.size_of_system()
  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
  A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

  iters = 0
  def nonlocal_iterate(vec_x):
    global iters
    iters += 1

  vectorSolution, num_iter = sp.linalg.cg(A, vectorRHS, rtol=1e-6, callback=nonlocal_iterate)
  if num_iter != 0:
    raise RuntimeError("Linear solver did not converge!")
  error = HDG_wrapper.errors(vectorSolution)[0]

  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_timoshenko_convergecne_cross.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Iteration = " + str(iteration) + 
          ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "timo_cross-" + str(poly_degree) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,7):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
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

