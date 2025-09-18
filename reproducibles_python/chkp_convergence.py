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
  
  goal_time = 1./1000.
  time_steps  = 1
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
    + str(3*poly_degree) + ",ChkpParametersZero,double>"
  const.cython_replacements = ["string", "string", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/chkp.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  lsol_constr = get_loc_constr(delta_time)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/square.geo", lsol_constr = get_loc_constr(delta_time) )
  HDG_wrapper.refine(iteration)

  vectorSolution = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  sol = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
  print(sol)
  print(HDG_wrapper.residual_flux(sol, goal_time))
  sol_res = (np.linalg.norm(HDG_wrapper.residual_flux(sol, goal_time)))**2/ len(vectorSolution)
  
  time = 0.
  #fun=lambda x: HDG_wrapper.residual_flux(x, time)
  def fun(x): 
    helper = (np.linalg.norm(HDG_wrapper.residual_flux(x, time)))**2/ len(vectorSolution)
    if helper < sol_res:
      sol = x.copy()
    return helper

  for time_step in range(time_steps):
    time += delta_time    
    #opt_obj = sp_opt.root(fun, vectorSolution, tol=1e-9, method='hybr', options ={"xtol": 1e-3})
    #opt_obj = sp_opt.root(fun, vectorSolution, tol=1e-6, method='krylov')
    opt_obj = sp_opt.minimize(fun, vectorSolution, method='BFGS', options={'xrtol': 1e-8, 'gtol': 1e-5})
    print(opt_obj.message)
    if not opt_obj.success:
      print(opt_obj.message)
      #raise RuntimeError("All linear solvers did not converge!")

    vectorSolution = opt_obj.x
    HDG_wrapper.set_data(vectorSolution, time)
    error = HDG_wrapper.errors(vectorSolution, time)[0]
    HDG_wrapper.set_data(sol, time)
    error2 = HDG_wrapper.errors(sol, time)[0]
    print( "Time: ", time, "\tError: ", error, "\t", error2 )
    
  error = HDG_wrapper.errors(vectorSolution, time)[0]
  print( "Iteration: ", iteration, " Error: ", error )
  
  #opt_obj = sp_opt.root(fun, vectorSolution, tol=1e-6, method='krylov')
  #opt_obj = sp_opt.minimize(fun, vectorSolution)
  #print(opt_obj.message)
  #vectorSolution = opt_obj.x
  
  #HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(2) + "-" + str(iteration) )
  #HDG_wrapper.plot_option( "printFileNumber" , "false" )
  #HDG_wrapper.plot_option( "scale" , "0.95" )
  #HDG_wrapper.plot_solution(vectorSolution, time + delta_time)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(0,1):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for iteration in [2, 4, 8, 16]:
      try:
        diffusion_test(poly_degree, iteration, debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

