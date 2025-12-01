from __future__ import print_function

import numpy as np
import scipy.optimize as sp_opt

import scipy.sparse as sp

from datetime import datetime

import os, sys

def get_loc_constr(h, t):
  return [t, 3., 3., 1., 1., 1., 1., 1., -3., -1., 3. + 1.]
  #Order: delta_t, tau+pu, tau-pu, tau-pv, tau+zu, tau-zu, tau-zv, tau+vu, tau_uqq, tau_yvu, tau_f


# --------------------------------------------------------------------------------------------------
# Function diffusion_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, iteration, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)
  os.system("mkdir -p output")
  
  h = 1. / iteration
  start_time  = 0.
  goal_time   = 4.
  time_steps  = 400

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
    + str(3*poly_degree) + ",ChkpParametersPeakon,double>"
  const.cython_replacements = ["string", "string", \
    "double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/chkp.hxx"]
  const.debug_mode      = debug_mode

  PyDP = HyperHDG.include(const)
  lsol_constr = get_loc_constr(h, delta_time)
  HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/lsq.geo", lsol_constr = get_loc_constr(h, delta_time) )
  HDG_wrapper.refine(iteration)
  
  def ttf_mat(x, time):
    col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(x, time)
    A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(len(x),len(x)))
    return A
  
  def reduce_shape(M):
    M.eliminate_zeros()
    M = M[M.getnnz(1) > 0]
    mask = M.getnnz(0) > 0
    M = M[:, mask]
    return M, mask

  def prolong(x, mask):
    r = np.zeros(mask.shape)
    r[mask] = x
    return r

  def newton(x, A, M, mask, time, tol=1e-8):
    rhs = np.array(HDG_wrapper.residual_flux(x, time))
    ra = np.linalg.norm(rhs)
    stepsize = 1.
    i = 0
    while ra > tol and i < 100:
      #print("Matrix assembliert")
      # A, mask = reduce_shape(A)
      # step = sp.linalg.gmres(A, rhs[mask], atol=1e-10, rtol=1e-2 * ra, M=sp.diags_array(1./A.diagonal()))[0]
      # step = prolong(step, mask)
      # print(datetime.now(), "Start solve")
      # sys.stdout.flush()
      step, info = sp.linalg.gmres(A, rhs[mask], atol=1e-10, rtol=1e-10, M=M)
      # print(info)
      step = prolong(step, mask)
      # print(datetime.now(), "End solve")
      # sys.stdout.flush()

      x   -= stepsize * step
      rhs  = np.array(HDG_wrapper.residual_flux(x, time))
      ra   = np.linalg.norm(rhs)
      i   += 1
      # print(datetime.now(),  i, ra)
      # sys.stdout.flush()
    return x

  time = start_time
  vectorSolution = np.array(HDG_wrapper.make_initial(HDG_wrapper.zero_vector(), time))

  for time_step in range(time_steps):
    time += delta_time
    if time_step % 10 == 0:
      # print(datetime.now(), "Start matrix")
      # sys.stdout.flush()
      A = ttf_mat(vectorSolution, time)
      A, mask = reduce_shape(A)
      # print(datetime.now(), "End matrix")
      # sys.stdout.flush()
      A_iLU = sp.linalg.spilu(A)
      M = sp.linalg.LinearOperator((np.sum(mask),np.sum(mask)), A_iLU.solve)
      # print(datetime.now(), "End preconditioner")
    x = newton(vectorSolution, A, M, mask, time)
    time = round(time, 8)
    
    res = np.linalg.norm(HDG_wrapper.residual_flux(vectorSolution, time))
    if round(time_steps * time / goal_time) % 1 == 0:
      HDG_wrapper.plot_option( "fileName" , "peakon" + str(poly_degree) + "-" + str(iteration) + "-" + str(time) )
      HDG_wrapper.plot_option( "printFileNumber" , "false" )
      HDG_wrapper.plot_option( "scale" , "1.0" )
      HDG_wrapper.plot_solution(vectorSolution, time)
  
    HDG_wrapper.set_data(vectorSolution, time)

    errors = HDG_wrapper.errors(vectorSolution, time)
    u_error = errors[0]
    q_error = errors[1]
    if round(time_steps * time / goal_time) % 1 == 0 or time == start_time + delta_time:
      print(datetime.now(), f'Time: {time:.6f}    Errors: {u_error:.2e} in u, {q_error:.2e} in q    Residual: {res}')
      sys.stdout.flush()
    
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in [2]:
    print("\nPolynomial degree is set to be ", poly_degree, "\n")
    for iteration in [128]:
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

