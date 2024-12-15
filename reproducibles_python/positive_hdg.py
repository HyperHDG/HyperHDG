from __future__ import print_function

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(debug_mode=False):
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
  const.topology        = "Cubic<2,2>"
  const.geometry        = "UnitCubes<2,2,double>"
  const.node_descriptor = "Cubic<2,2>"
  const.local_solver    = "DiffusionSpaces<2,0,1,0,4,TestParametersPos,double>"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( [1,1],[1,1], 1. )

  vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

  system_size = HDG_wrapper.size_of_system()
  # A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
  A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

  [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, rtol=1e-13)
  if num_iter != 0:
    print("CG solver failed with a total number of ", num_iter, "iterations.")
    [vectorSolution, num_iter] = sp_lin_alg.gmres(A, vectorRHS, rtol=1e-13)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations.")
      raise RuntimeError("Linear solvers did not converge!")

  print(vectorSolution)
  
  HDG_wrapper.plot_option( "fileName" , "positive_hdg" )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  # HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  diffusion_test(debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
