from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime
import argparse
import os, sys

parser = argparse.ArgumentParser(prog="heat")
parser.add_argument("-i", type=int, default=1)
parser.add_argument("-p", type=int, default=3)
parser.add_argument("-n", type=int, default=10)
parser.add_argument("-d", type=int, default=1)
parser.add_argument("-t", type=float, default=1.)
parser.add_argument("-o")

args = parser.parse_args()
iteration = args.i
poly_degree = args.p
dimension = args.d
debug_mode = False

theta       = args.t
time_steps  = 2 ** args.n
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
    [vectorSolution, num_iter] = sp_lin_alg.gmres(A,vectorRHS,rtol=1e-13)
    if num_iter != 0:
      print("GMRES also failed with a total number of ", num_iter, "iterations.")
      [vectorSolution, num_iter] = sp_lin_alg.bicgstab(A,vectorRHS,rtol=1e-13)
      if num_iter != 0:
        print("BiCGStab also failed with a total number of ", num_iter, "iterations.")
        raise RuntimeError("All linear solvers did not converge!")

  HDG_wrapper.set_data(vectorSolution, (time_step+1) * delta_time)
  
error = HDG_wrapper.errors(vectorSolution, 1.)[0]

print("# heat.py")
print(f"theta: {theta}")
print(f"dimension: {dimension}")
print(f"degree: {poly_degree}")
print(f"iteration: {iteration}")
print(f"timesteps: {time_steps}")
print(f"error: {error:.5e}")

if args.o:
    HDG_wrapper.plot_option( "fileName" , args.o)
    HDG_wrapper.plot_option( "printFileNumber" , "false" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, 1.)
    print(f"output: output/{args.o}.vtu")
