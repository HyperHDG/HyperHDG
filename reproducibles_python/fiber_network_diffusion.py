from __future__ import print_function

import numpy as np
import scipy.sparse as sp

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, debug_mode=False):
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
  const.topology        = "File<1,2,std::vector,Point<2,double> >"
  const.geometry        = "File<1,2,std::vector,Point<2,double> >"
  const.node_descriptor = "File<1,2,std::vector,Point<2,double> >"
  const.local_solver    = "Diffusion<1," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",TestParametersSinEllipt,double>"
  const.cython_replacements = ["string", "string"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/fiber_network_1000.geo" )

  vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

  system_size = HDG_wrapper.size_of_system()
  # A = sp.linalg.LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
  A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))


  points = np.loadtxt("domains/fiber_network_1000_points.txt")
  helper = HyperHDG.fiber_network.precond(points, [2**3, 2**3])
  def precond_mult( vec_x ):
    return helper.precond(A, vec_x)
  B = sp.linalg.LinearOperator( (system_size,system_size), matvec= precond_mult )


  iters = 0
  def nonlocal_iterate(vec_x):
    nonlocal iters
    iters += 1
    print(iters, "\t", np.linalg.norm(A.dot(vec_x) - vectorRHS) / np.linalg.norm(vectorRHS),
          "\t", .5 * vec_x.dot(A.dot(vec_x)) - vec_x.dot(vectorRHS), " \t", datetime.now())

  [vectorSolution, n_iter] = sp.linalg.cg(A, vectorRHS, rtol=1e-10, callback=nonlocal_iterate, M=B)
  print(iters)
  if n_iter != 0:
    print("CG solver failed with a total number of ", n_iter, "iterations.")
    [vectorSolution, n_iter] = sp.linalg.gmres(A, vectorRHS)
    if n_iter != 0:
      print("GMRES also failed with a total number of ", n_iter, "iterations.")
      raise RuntimeError("Linear solvers did not converge!")

  error = HDG_wrapper.errors(vectorSolution)[0]
  print(" Error: ", error)
  f = open("output/diffusion_convergence_elliptic.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Error = " + str(error) + ".\n")
  f.close()
  
  HDG_wrapper.plot_option( "fileName" , "fiber_network_1000-" + str(poly_degree) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution)

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in [1, 6]:
    print("\n Polynomial degree is set to be ", poly_degree, "\n")
    try:
      diffusion_test(poly_degree, debug_mode)
    except RuntimeError as error:
      print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
