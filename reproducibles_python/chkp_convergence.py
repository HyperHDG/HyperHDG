from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

import scipy.sparse as sp

from datetime import datetime

import os, sys

def get_loc_constr(h, t):
  return [t, 3. + 1. / h, 3. + 1. / h, 1., 1., 1., 1., 1., -3., -2]


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  h = 1. / iteration
  goal_time = .5 
  time_steps  = 16 * iteration

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
  lsol_constr = get_loc_constr(h, delta_time)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/square.geo", lsol_constr = get_loc_constr(h, delta_time) )
  HDG_wrapper.refine(iteration)
  
  def ttf_mat(x, time):
    col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
    A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
    return A

  def rf(x):
    return np.array(HDG_wrapper.residual_flux(x, time))

  def newton(x, time, tol=1e-8):
    ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time))
    stepsize = 1.
    i = 0
    while ra > tol and i < 100:
      A = ttf_mat(x, time)
      step = sp.linalg.gmres(A, HDG_wrapper.residual_flux(x, time), atol=1e-10, rtol=1e-10)[0]
      x -= stepsize * step
      ra = np.linalg.norm(HDG_wrapper.residual_flux(x, time))
      i += 1
      print(i, ra)
    return x

  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector()))
  time = 0.

  for time_step in range(time_steps):
    time += delta_time
    #newton(vectorSolution, time)
    # opt_obj = sp_opt.root(rf, vectorSolution, jac=lambda x:ttf_mat(x, time).todense(), method='hybr')
    # vectorSolution = opt_obj.x
    x = newton(vectorSolution, time)
    #time = round(time, 8)
    newton(vectorSolution, time)
    #opt_obj = sp_opt.root(rf, vectorSolution, jac=lambda x:ttf_mat(x, time).todense(), method='hybr')
    #vectorSolution = opt_obj.x
    
    res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
    #print(vectorSolution)
    HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(poly_degree) + "-" + str(iteration) + "-" + str(time) )
    HDG_wrapper.plot_option( "printFileNumber" , "false" )
    #HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, time)
  
    HDG_wrapper.set_data(vectorSolution, time)

    u_error = HDG_wrapper.errors(vectorSolution, time)[0]
    q_error = HDG_wrapper.errors(vectorSolution, time)[1]
    print(f'{f'Time: {time:.6f}':20}Errors: {u_error:.2e} in u, {q_error:.2e} in q\tResidual: {res}')
    sys.stdout.flush()
    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for iteration in [4, 8, 16, 32]:
    print("\n\n Grid size is set to be ", iteration)
    for poly_degree in [2, 3]:
      print("\nPolynomial degree is set to be ", poly_degree, "\n")
      try:
        diffusion_test(poly_degree, iteration, debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

