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
def analytic_test():
  poly_degree = 1
  iteration = 2
  debug_mode = False

  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
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
    + str(3*poly_degree) + ",ChkpParametersLinear,double>"
  const.cython_replacements = ["string", "string", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/chkp.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  lsol_constr = get_loc_constr(1.)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/unitsquare.geo", lsol_constr = get_loc_constr(1.) )
  HDG_wrapper.refine(iteration)
      

  vsr = [-1.,         0.,         1.,         0.,         0.,         0.,\
 -1.,         0.,         1.,         0.,         0.,         0.,\
 0.,         0.,         1.,         0.,         0.,         0.,\
 0.,         0.,         1.,         0.,         0.,         0.,\
 -0.75,      0.14433757, 0.,        0.,        0.,        0.,\
 -0.25,      0.14433757, 0.,        0.,        0.,        0.,\
 -0.75,      0.14433757, 0.,        0.,        0.,        0.,\
 -0.25,      0.14433757, 0.,        0.,        0.,        0.,\
 -0.5,       0.,        1.,        0.,        0.,        0.,\
 -0.5,       0.,        1.,        0.,        0.,        0.,\
 -0.75,      0.14433757, 0.,        0.,        0.,        0.,\
 -0.25,      0.14433757, 0.,        0.,        0.,        0.        ]  
  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector())) 
  vectorSolution = vsr

  time = 1

  res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))

  HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(poly_degree) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time)
  
  HDG_wrapper.set_data(vectorSolution, time)
  error = HDG_wrapper.errors(vectorSolution, time)[0]
  print( "Time: ", time, "\tError: ", error, "\t", " Residual: ", res)

    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  try:
    analytic_test()
  except RuntimeError as error:
    print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

