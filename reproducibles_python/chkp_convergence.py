from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

import scipy.sparse as sp

from datetime import datetime

import os, sys

def get_loc_constr(h, t):
  return [t,       -1.,     -1.,     1.,      -2,       4.   ]
  #Order: delta_t, tau+zpu, tau-zpu, tau-zpv, tau_uqq,  tau_f


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  h = 1. / iteration
  goal_time = .1
  time_steps  = 100

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
  
  # def ttf_mat(x, time):
  #   col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
  #   A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
  #   At = sp.csc_matrix((vals, (col_ind, row_ind)), shape=(len(x),len(x)))
  #   return At @ A, At
  
  def ttf_mat(x, time):
    col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
    A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
    return A
  
  # def reduce_shape(M):
  #   M.eliminate_zeros()
  #   mask_r = M.getnnz(1) > 0
  #   M = M[mask_r]
  #   mask_c = M.getnnz(0) > 0
  #   M = M[:, mask_c]
  #   return M, mask_r, mask_c


  # Alternative version that assumes square matrix and removes rows/columns at same indices
  def remove_zero_rows_and_columns(csc_matrix):
    # Find indices of non-zero rows
    col_sums = csc_matrix.sum(axis=0).A1
    row_sums = csc_matrix.sum(axis=1).A1

    non_zero_cols = np.where((col_sums != 0))[0]
    non_zero_rows = np.where((row_sums != 0))[0]
        
    csc_matrix = csc_matrix[non_zero_rows, :]
    csc_matrix = csc_matrix[:, non_zero_cols]
    
    return csc_matrix, non_zero_cols, non_zero_rows



  # def prolong(x, mask_c):
  #   r = np.zeros(mask_c.shape)
  #   r[mask_c] = x
  #   return r

  def prolong(x, keep_rows, full_length):
    r = np.zeros(full_length,)
    r[keep_rows] = x
    return r

  # def rf(x):
  #   return np.array(HDG_wrapper.residual_flux(x, time))

  # def newton(An, At, mask_c, x, time, tol=1e-8):
  #   rhs = np.array(HDG_wrapper.residual_flux(x, time))[mask_c]
  #   ra = np.linalg.norm(rhs)
  #   stepsize = 1.
  #   i = 0
  #   while ra > tol and i < 100:
  #     step = sp.linalg.gmres(An, At @ rhs, atol=1e-10, rtol=1e-10)[0]
  #     step = prolong(step, mask_c)
  #     # print(datetime.now(), "End solve")
  #     # sys.stdout.flush()

  #     x   -= stepsize * step
  #     rhs  = np.array(HDG_wrapper.residual_flux(x, time))
  #     ra   = np.linalg.norm(rhs)
  #     rhs = rhs[mask_c]
  #     i   += 1
  #     print(datetime.now(),  i, ra)
  #     # sys.stdout.flush()
  #   return x

  def newton(A, M, keep_cols, keep_rows, x, time, tol=1e-8):
    rhs  = np.array(HDG_wrapper.residual_flux(x, time))
    rhs_len = len(rhs)
    rhs  = rhs[keep_rows]
    for _ in range(10):
      step, _ = sp.linalg.gmres(A, rhs, M=M, atol=1e-10, rtol=1e-10)
      x   -= prolong(step, keep_cols, rhs_len)
      rhs = np.array(HDG_wrapper.residual_flux(x, time))
      rhs = rhs[keep_rows]
      if np.linalg.norm(rhs) < tol:  return x
    print("Newton failed!")
    return x

  # vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector()))
  # An, At = ttf_mat(vectorSolution, delta_time)
  # #An, mask_r, mask_c = reduce_shape(An)
  # #At = At[mask_r]
  # #At = At[:, mask_c]
  # #M = sp.linalg.LinearOperator(An.shape, sp.linalg.spilu(An).solve)
  # mask_c = np.ones(An.shape[0], dtype='bool')
  # time = 0.

  time = 0.
  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector()))

  # for time_step in range(time_steps):
  #   time += delta_time
  #   if time_step % 10 == 1:
  #     An, At = ttf_mat(vectorSolution, time)
  #     #An, mask_r, mask_c = reduce_shape(An)
  #     #At = At[mask_r]
  #     #At = At[:, mask_c]
  #     #M = sp.linalg.LinearOperator(An.shape, sp.linalg.spilu(An).solve)
  #   vectorSolution = newton(An, At, mask_c, vectorSolution, time)
  #   time = round(time, 8)

  for time_step in range(time_steps):
    time += delta_time
    if time_step % 10 == 0:
      A = ttf_mat(vectorSolution, time)
      A, keep_cols, keep_rows = remove_zero_rows_and_columns(A)
      assert len(keep_cols) == len(keep_rows), "Error in removing zero rows and columns!"
      sA_iLU = sp.linalg.spilu(A)
      M = sp.linalg.LinearOperator((len(keep_rows),len(keep_rows)), sA_iLU.solve)

    vectorSolution = newton(A, M, keep_cols, keep_rows, vectorSolution, time)
    time = round(time, 8)
    
    res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
    HDG_wrapper.plot_option( "fileName" , "chkp_conv" + str(poly_degree) + "-" + str(iteration) + "-" + str(time) )
    HDG_wrapper.plot_option( "printFileNumber" , "false" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, time)
  
    HDG_wrapper.set_data(vectorSolution, time)

    errors = HDG_wrapper.errors(vectorSolution, time)
    u_error = errors[0]
    q_error = errors[1]
    # if round(time_steps * time / goal_time) % 1 == 0 or time == delta_time:
    if (time_step+1) % 10 == 0:
      print(datetime.now(), f'Time: {time:.6f}    Errors: {u_error:.2e} in u, {q_error:.2e} in q    Residual: {res}')
      sys.stdout.flush()
    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in [1, 2, 3]:
    print("\nPolynomial degree is set to be ", poly_degree, "\n")
    for iteration in [2, 4, 8, 16]:
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

