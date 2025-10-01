from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys
import argparse
import prin2
import logging

# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration, debug_mode=False):
  logger = logging.getLogger("diffusion_parabolic")

  os.system("mkdir -p output")
  
  theta       = 1.
  time_steps  = 10 ** 4
  delta_time  = 1 / time_steps
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Parabolic"
  const.topology        = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
  const.geometry        = "UnitCube<" + str(dimension) + "," + str(dimension) + ",double>"
  const.node_descriptor = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
  const.local_solver    = "DiffusionParab<" + str(dimension) + "," + str(poly_degree) + "," \
    + str(2*poly_degree) + ",TestParametersSinParab,double>"
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  HDG_wrapper = PyDP( [2 ** iteration] * dimension, lsol_constr= [1.,theta,delta_time] )

  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  
  system_size = HDG_wrapper.size_of_system()
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

  for time_step in range(time_steps):
    
    vectorRHS = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
                 (time_step+1) * delta_time), -1.)
    
    [vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, rtol=1e-13)
    if num_iter != 0:
      print("CG failed with a total number of ", num_iter, " iterations in time step ", time_step, \
            ". Trying GMRES!")
      [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,tol=1e-13)
      if num_iter != 0:
        print("GMRES also failed with a total number of ", num_iter, "iterations.")
        [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,tol=1e-13)
        if num_iter != 0:
          print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
          raise RuntimeError("All linear solvers did not converge!")

    HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)
    
  error = HDG_wrapper.errors(vectorSolution, 1.)[0]
  logger.info(f"{iteration=}, {error=}")
  
  HDG_wrapper.plot_option( "fileName" , "diff_conv_parab" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, 1.)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
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
  parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
  parser.add_argument("-d", "--dimension", help="dimension of the problem", default=1, type=int)
  parser.add_argument("-p", "--degree",    help="polynomial degree of approximation", default=1, type=int)
  parser.add_argument("-i", "--iteration", help="iteration", default=1, type=int)
  parser.add_argument("--debug", help="toggle debug mode", action="store_true")
  parser.add_argument("--log-level", help="set the log level")
  args = parser.parse_args()

  logging.setLoggerClass(prin2.Logger)
  log_levels = {
      'debug': logging.DEBUG,
      'info': logging.INFO,
      'warning': logging.WARNING,
      'error': logging.ERROR,
      'critical': logging.CRITICAL
  }
  logger = logging.getLogger("diffusion_parabolic")
  logger.setLevel(level=log_levels.get(args.log_level, logging.INFO))
  logger.log_args(args)

  diffusion_test(args.degree, args.dimension, args.iteration, args.debug)
