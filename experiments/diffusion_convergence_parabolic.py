from __future__ import print_function

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator

from datetime import datetime, timedelta

import os, sys
import argparse
import prin2
import logging

parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("-d", "--dimension", help="dimension of the problem", default=1, type=int)
parser.add_argument("-p", "--degree",    help="polynomial degree of approximation", default=1, type=int)
parser.add_argument("-i", "--iteration", help="iteration", default=1, type=int)
parser.add_argument("--debug", help="toggle debug mode", action="store_true")
parser.add_argument("--log-level", help="set the log level")
parser.add_argument("--rtol", help="rtol", type=float, default=1e-10)
parser.add_argument("--direct", help="direct", action="store_true")
parser.add_argument("--mat", help="mat", action="store_true")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-n", "--time-steps", help="output file", type=int, default=10**2)
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

poly_degree = args.degree
dimension = args.dimension
iteration = args.iteration
debug_mode = args.debug

logger = logging.getLogger("diffusion_parabolic")

theta       = 1.
time_steps  = args.time_steps
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
  + str(2*poly_degree) + ",TestHeat,double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", \
  "double", "vector[double]"]
const.include_files   = ["reproducibles_python/parameters/diffusion.hxx", "experiments/parameters.hxx"]
const.debug_mode      = debug_mode

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP( [2 ** iteration] * dimension, lsol_constr= [1.,theta,delta_time] )


system_size = HDG_wrapper.size_of_system()
logger.info(f"{system_size=}")

if args.direct or args.mat:
  logger.info("assembling matrix...")
  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
  A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
  A.eliminate_zeros()
  print(A)
  if args.direct:
    splu = sp.linalg.splu(A)
else:
  A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )
  
iters = 0
start = datetime.now()
def count_iter(x):
  global iters
  iters +=1

logger.info("timestepping...")

initial = HDG_wrapper.make_initial(HDG_wrapper.zero_vector())
if args.output:
    HDG_wrapper.plot_option( "fileName" , args.output)
    HDG_wrapper.plot_option( "printFileNumber" , "true" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(initial, 0.)

for time_step in range(time_steps):
  rhs = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector(), \
               (time_step+1) * delta_time), -1.)
  
  [vectorSolution, num_iter] = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=count_iter)

  if num_iter != 0:
    logger.error(f"no convergence in {num_iter} iterations")
    break

  if args.output:
    HDG_wrapper.set_data(vectorSolution, (time_step+1)*delta_time)
    HDG_wrapper.plot_option( "fileName" , args.output)
    HDG_wrapper.plot_option( "printFileNumber" , "true" )
    HDG_wrapper.plot_option( "scale" , "0.95" )
    HDG_wrapper.plot_solution(vectorSolution, (time_step+1)*delta_time)
  
error = HDG_wrapper.errors(vectorSolution, 1.)[0]
logger.info(f"{iteration=}, {error=}, avg num iters={iters/time_steps}")
