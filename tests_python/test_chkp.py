from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

import scipy.sparse as sp

from datetime import datetime

import os, sys

def get_loc_constr(t):
  return [t, 3., 4., 1., 1., 1., 1., 1., -3., -2]


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test():
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

  def newton(x, time, tol=1e-8):
    ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time))
    stepsize = 1.
    i = 0
    while ra > tol and i < 100:
      col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
      A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
      step = sp.linalg.lsqr(A, HDG_wrapper.residual_flux(x, time))[0]
      x -= stepsize * step
      ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time))
      i += 1
      print(i, ra)
    return x
      

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
  print(vectorSolution)
  #print(sol)
  #vectorSolution += 0. * np.random.uniform(-1., 1., vectorSolution.shape)
  #print(HDG_wrapper.residual_flux(vectorSolution, goal_time))

  time = 1
#  fun=lambda x: HDG_wrapper.residual_flux(x, time)
#  def fun(x): 
#    helper = (np.linalg.norm(HDG_wrapper.residual_flux(x, time)))**2/ len(vectorSolution)
#    return helper
  def fun(x):
    helper = HDG_wrapper.residual_flux(x, time)
#    print(np.linalg.norm(helper))
    return helper

  newton(vectorSolution, time)

  res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
  print(HDG_wrapper.residual_flux(vectorSolution, time))
    #opt_obj = sp_opt.root(fun, vectorSolution, tol=1e-9, method='hybr', options ={"xtol": 1e-3})
    #opt_obj = sp_opt.root(fun, vectorSolution, tol=1e-9, method='krylov')
    #opt_obj = sp_opt.minimize(fun, vectorSolution, method='BFGS', options={'xrtol': 1e-8, 'gtol': 1e-5})
    #print(opt_obj.message, opt_obj.nit, np.linalg.norm(opt_obj.x))
    #if not opt_obj.success:
    #  print(opt_obj.message)
    #  #raise RuntimeError("All linear solvers did not converge!")

    #vectorSolution = opt_obj.x

  deriv_err = 0.
  for i in []:
    vec_dir = np.array(HDG_wrapper.zero_vector())
    vec_dir[i] = 1.
    deriv_diff = HDG_wrapper.trace_to_flux(vectorSolution, vec_dir, 1.)
    #print(deriv_diff)
    finite_diff = (np.array(HDG_wrapper.residual_flux(vectorSolution + 1e-4 * vec_dir, 1.)) - np.array(HDG_wrapper.residual_flux(vectorSolution, 1.))) / 1e-4
    #print(finite_diff)
    deriv_diff -= finite_diff
    #print(deriv_diff)
    deriv_err = max(np.max(np.abs(deriv_diff)), deriv_err)
    #print(i, deriv_err)

  HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(poly_degree) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution, time)
  
  HDG_wrapper.set_data(vectorSolution, time)
  error = HDG_wrapper.errors(vectorSolution, time)[0]
  print( "Time: ", time, "\tError: ", error, "\t", " Residual: ", res)
  print(vectorSolution - vsr)

    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  try:
    diffusion_test()
  except RuntimeError as error:
    print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

