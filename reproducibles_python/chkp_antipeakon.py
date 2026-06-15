from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

import scipy.sparse as sp
import scipy.optimize as sp_opt

from datetime import datetime

import os, sys

import warnings
warnings.filterwarnings("error", category=RuntimeWarning)

def get_loc_constr(h, t):
  return [t,       -1.,     -1.,     1.,      -0.25,    4.   ]
  #Order: delta_t, tau+zpu, tau-zpu, tau-zpv, tau_uqq,  tau_f


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, iteration, debug_mode=False):
  begin_time = datetime.now()
  print("Starting time is", begin_time)
  os.system("mkdir -p output")
  
  h = 1. / iteration
  start_time  = 0.00000
  goal_time   = 0.00001
  time_steps  = 1

  delta_time  = (goal_time - start_time) / time_steps
  
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
    + str(3*poly_degree) + ",ChkpParametersAntipeakon,double>"
  const.cython_replacements = ["string", "string", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/chkp.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  lsol_constr = get_loc_constr(h, delta_time)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/lsq2.geo", lsol_constr = get_loc_constr(h, delta_time) )
  HDG_wrapper.refine(iteration)
  
  def ttf_mat(x, time):
    col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
    A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
    return A

  def remove_zero_rows_and_columns(csc_matrix):
    # Find indices of non-zero rows
    col_sums = csc_matrix.sum(axis=0).A1
    row_sums = csc_matrix.sum(axis=1).A1

    non_zero_cols = np.where((col_sums != 0))[0]
    non_zero_rows = np.where((row_sums != 0))[0]
        
    csc_matrix = csc_matrix[non_zero_rows, :]
    csc_matrix = csc_matrix[:, non_zero_cols]
    
    return csc_matrix, non_zero_cols, non_zero_rows


  def prolong(x, keep_rows, full_length):
    r = np.zeros(full_length,)
    r[keep_rows] = x
    return r

  def newton(A, M, keep_cols, keep_rows, x, time, tol=1e-10):
    rhs  = np.array(HDG_wrapper.residual_flux(x, time))
    rhs_len = len(rhs)
    rhs  = rhs[keep_rows]
    for _ in range(25):
      step, _ = sp.linalg.gmres(A, rhs, M=M, atol=1e-10, rtol=1e-10)
      x   -= prolong(step, keep_cols, rhs_len)
      rhs = np.array(HDG_wrapper.residual_flux(x, time))
      rhs = rhs[keep_rows]
      if np.linalg.norm(rhs) < tol:  return x
    raise ValueError("Frozen Newton failed")

  def fallback_newton(x, time, tol=1e-8, alpha=0.1, beta=0.5):
    A = ttf_mat(x, time)
    A, keep_cols, keep_rows = remove_zero_rows_and_columns(A)
    assert len(keep_cols) == len(keep_rows), "Error in removing zero rows and columns!"
    sA_iLU = sp.linalg.spilu(A)
    M = sp.linalg.LinearOperator((len(keep_rows),len(keep_rows)), sA_iLU.solve)
    rhs  = np.array(HDG_wrapper.residual_flux(x, time))
    rhs_len = len(rhs)
    rhs  = rhs[keep_rows]
    res_old = np.linalg.norm(rhs)
    for _ in range(50):
      step, info = sp.linalg.gmres(A, rhs, M=M, atol=1e-10, rtol=1e-10)
      print(res_old, info, np.linalg.norm(A@step - rhs))
      x_cand = x - prolong(step, keep_cols, rhs_len)
      stepsize = 1.
      rhs = np.array(HDG_wrapper.residual_flux(x_cand, time))
      rhs = rhs[keep_rows]
      res = np.linalg.norm(rhs)
      while stepsize > 1e-16 and res > (1. - alpha * stepsize) * res_old:
        stepsize *= beta
        x_cand = x - stepsize * prolong(step, keep_cols, rhs_len)
        rhs = np.array(HDG_wrapper.residual_flux(x_cand, time))
        rhs = rhs[keep_rows]
        res = np.linalg.norm(rhs)
      if stepsize <= 1e-16:
        raise ValueError("Stepsize too small")
      x = x_cand
      res_old = res
      if res_old < tol:
        return A, M, keep_cols, keep_rows, x
      A = ttf_mat(x, time)
      A, keep_cols, keep_rows = remove_zero_rows_and_columns(A)
      assert len(keep_cols) == len(keep_rows), "Error in removing zero rows and columns!"
      sA_iLU = sp.linalg.spilu(A)
      M = sp.linalg.LinearOperator((len(keep_rows),len(keep_rows)), sA_iLU.solve)
    raise ValueError("Fallback Newton failed")
    return A, M, keep_cols, keep_rows, x

  def fallback_newton_step(x, time):
    A = ttf_mat(x, time)
    A, keep_cols, keep_rows = remove_zero_rows_and_columns(A)
    assert len(keep_cols) == len(keep_rows), "Error in removing zero rows and columns!"
    sA_iLU = sp.linalg.spilu(A)
    M = sp.linalg.LinearOperator((len(keep_rows),len(keep_rows)), sA_iLU.solve)
    rhs  = np.array(HDG_wrapper.residual_flux(x, time))
    rhs_len = len(rhs)
    rhs  = rhs[keep_rows]
    step, _ = sp.linalg.gmres(A, rhs, M=M, atol=1e-10, rtol=1e-10)
    x   -= prolong(step, keep_cols, rhs_len)
    return A, M, keep_cols, keep_rows, x

  time = start_time
  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector(), time))

  for time_step in range(time_steps):
    time += delta_time

    if time_step == 0:
      A, M, keep_cols, keep_rows, vectorSolution = fallback_newton_step(vectorSolution, time)
      print("First additional step")

    vs = np.copy(vectorSolution)
    try:
      vectorSolution = newton(A, M, keep_cols, keep_rows, vectorSolution, time)
    except (ValueError, RuntimeWarning, np.linalg.LinAlgError) as e:
      print(f'Compute new matrix at time {time:.6f}: ', e)
      vectorSolution = np.copy(vs)
      A, M, keep_cols, keep_rows, vectorSolution = fallback_newton_step(vectorSolution, time)
      try:
        vectorSolution = newton(A, M, keep_cols, keep_rows, vectorSolution, time)
      except (ValueError, RuntimeWarning, np.linalg.LinAlgError) as e:
        print("Try fallback Newton instead: ", e)
        vectorSolution = np.copy(vs)
        A, M, keep_cols, keep_rows, vectorSolution = fallback_newton(vectorSolution, time)
        '''
        try:
          A, M, keep_cols, keep_rows, vectorSolution = fallback_newton(vectorSolution, time)
        except ValueError as e:
          print("Fallback Newton failed, try scipy root instead")
          vectorSolution = np.copy(vs)
          fun = lambda x: np.array(HDG_wrapper.residual_flux(x, time))
          jac = lambda x: ttf_mat(x, time).todense()
          res = sp_opt.root(fun, vectorSolution, tol=1e-6, jac=jac, method='lm')
          vectorSolution = res.x
          print(f'Residual norm is {np.linalg.norm(res.fun)}')
          A = ttf_mat(vectorSolution, time)
          A, keep_cols, keep_rows = remove_zero_rows_and_columns(A)
          assert len(keep_cols) == len(keep_rows), "Error in removing zero rows and columns!"
          sA_iLU = sp.linalg.spilu(A)
          M = sp.linalg.LinearOperator((len(keep_rows),len(keep_rows)), sA_iLU.solve)
          '''



      
      
      
    if ((time_step - 1) % 200 == 0):
      A, M, keep_cols, keep_rows, vectorSolution = fallback_newton_step(vectorSolution, time)

    time = round(time, 8)
    
    res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
    if (time_step+1) % 1 == 0:
      HDG_wrapper.plot_option( "fileName" , "antipeakon_v" + str(poly_degree) + "-" + str(iteration) + "-" + str(time) )
      HDG_wrapper.plot_option( "printFileNumber" , "false" )
      HDG_wrapper.plot_option( "scale" , "1.0" )
      HDG_wrapper.plot_solution(vectorSolution, time)
  
    HDG_wrapper.set_data(vectorSolution, time)

    errors = HDG_wrapper.errors(vectorSolution, time)
    u_error = errors[0]
    q_error = errors[1]
    if ((time_step+1) % 100 == 0 and abs(time - 2.5) < 0.5) or (time_step + 1) % 100 == 0 or time_step == 0:
      print(datetime.now(), f'Time: {time:.6f}    Errors: {u_error:.2e} in u, {q_error:.2e} in q    Residual: {res}')
      sys.stdout.flush()
    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-begin_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in [2]:
    print("\nPolynomial degree is set to be ", poly_degree, "\n")
    for iteration in [32]:
      print("\n\n Grid size is set to be ", iteration)
      try:
        diffusion_test(poly_degree, iteration, debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")

