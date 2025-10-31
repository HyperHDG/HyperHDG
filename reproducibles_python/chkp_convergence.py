from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

from datetime import datetime

import os, sys

def get_loc_constr(t):
  return [t, 3., 4., 1., 1., 1., 1., 1., -3., -2]


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  goal_time = 1e-2
  time_steps  = 10
  delta_time  = goal_time / time_steps
  
  try:
    import HyperHDG
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
    import HyperHDG
  
  const                 = HyperHDG.config()
  const.global_loop     = "Nonlinear"
  const.topology        = "File<2,2>"
  const.geometry        = "File<2,2>"
  const.node_descriptor = "File<2,2>"
  const.local_solver    = "Chkp<" + str(2) + "," + str(poly_degree) + "," \
    + str(3*poly_degree) + ",ChkpParameters,double>"
  const.cython_replacements = ["string", "string", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/chkp.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  lsol_constr = get_loc_constr(delta_time)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/square.geo", lsol_constr = get_loc_constr(delta_time) )
  HDG_wrapper.refine(iteration)
  
  def jacobi(x, time):
    def partial(x, i, time):
      dv = np.array(HDG_wrapper.zero_vector())
      dv[i] = 1.
      return HDG_wrapper.trace_to_flux(x, dv, time)
    return np.array([partial(x, i, time) for i in range(len(x))]).transpose()

  def newton(x, time, tol=1e-10):
    ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time)) / len(x)
    stepsize = 1.
    i = 0
    while ra > tol and i < 100:
      step = np.linalg.lstsq(jacobi(x, time), HDG_wrapper.residual_flux(x, time))[0]
      x -= stepsize * step
      ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time)) / len(x)
      i += 1
      print(i, ra)
    return x

  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector()))
  time = 0.

  for time_step in range(time_steps):
    time += delta_time
    newton(vectorSolution, time)
    
    res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
    HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(poly_degree) + "-" + str(iteration) + "-" + str(time) )
    HDG_wrapper.plot_option( "printFileNumber" , "false" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, time + delta_time)
  
    HDG_wrapper.set_data(vectorSolution, time)
    error = HDG_wrapper.errors(vectorSolution, time)[0]
    print( "Time: ", time, "\tError: ", error, "\t", " Residual: ", res)

    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in [1]:
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for iteration in [1, 2, 4, 8]:
      try:
        diffusion_test(poly_degree, iteration, debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

