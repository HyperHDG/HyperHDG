#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
# The data necessary to run this file can be obtained using the HyperHDG.fiber_network.make_geo
# function on the files of Hauck, M., & Rupp, A. (2024). Fiber network models of paper. Zenodo.
# https://doi.org/10.5281/zenodo.12751486
# --------------------------------------------------------------------------------------------------

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import os, sys, argparse, logging, datetime, subprocess
import zstandard as zstd

import prin2, jprecond, geobin

######### SETUP argument parsing & logging

default_output_dir = "output"
now = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
default_output_name = f"{os.path.basename(__file__)}.{now}"
parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("network", help="full path to the network file")
parser.add_argument("domains", help="domains file")
parser.add_argument("-t", "--rtol",
  help="relative tolerance when to stop the CG iterator",
  type=float, default=1e-10
)
parser.add_argument("-d", "--debug", help="toggle debug mode", action="store_true")
parser.add_argument("-o", "--output", help="output name", default=default_output_name)
parser.add_argument("--output-dir", help="output dir", default=default_output_dir)
parser.add_argument("--maxiter", help="maximum number of cg iterations", type=int, default=100)
parser.add_argument("-m", "--modelproblem", help="the model problem to select", default="timo")
parser.add_argument("-n","--num-elements", help="the number of elements to use in the coarse finite element mesh", default=2**3, type=int)
parser.add_argument("--mat", help="store/load the lhs HDG matrix, depending on wether the path exists")

logging.setLoggerClass(prin2.Logger)
logger = logging.getLogger("fiber_network_elastic")

args = parser.parse_args()
logger.log_args(args)

# verify output path is writable
output_path = f"{args.output_dir}/{args.output}"
try:
  with open(output_path + ".vtu", "w") as file:
    file.write("test")
except Exception as e:
  print(f"error: cannot write to output path '{output_path}")
  sys.exit(1)

######## MAIN code

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG
  
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = args.debug
repeat = 1

if args.modelproblem == "timo":
  const.local_solver = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>"
  repeat = 6 # for every node we have 6 unknowns, namely displacement+rotation
elif args.modelproblem == "diff":
  const.local_solver    = "Diffusion<1,5,10,ConstantDiffusionParameters>"
  const.include_files   = ["experiments/parameters.hxx"]

logger.info("reading files")

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP(args.network)

network_points = geobin.read_network_points(args.network)

domains = None
if args.domains:
  domains = geobin.read_domains(args.domains)

logger.info("computing residual")

rhs = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

logger.info("assembling  A...")

system_size = HDG_wrapper.size_of_system()
if args.mat and os.path.isfile(args.mat):
  logger.info(f"  read matrix A from '{args.mat}'")
  A = sp.load_npz(args.mat)
else:
  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
  A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
  sp.save_npz(args.mat, A)
  logger.info(f"  wrote matrix A to '{args.mat}'")

logger.info("assembling  B...")

precond = jprecond.JPrecond(
  A,
  network_points,
  [args.num_elements, args.num_elements],
  repeat=repeat,
  domains=domains
)
B = sp.linalg.LinearOperator(
  (system_size,system_size),
  matvec=precond.matmul
)

iters = 0
start = datetime.datetime.now()
avg_time = datetime.timedelta(0)

def log_iter(x):
  global iters, start, avg_time
  iters += 1
  relErr = np.linalg.norm(A @ x - rhs) / np.linalg.norm(rhs)
  bilin = .5 * x.dot(A @ x) - x.dot(rhs)
  duration = datetime.datetime.now() - start
  start = datetime.datetime.now()
  avg_time += duration
  logger.info(f"{iters:>5} {relErr:>13.6e} {bilin:>13.6e} {duration}")

logger.info("starting cg")

vectorSolution, num_iter = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=log_iter, M=B, maxiter=args.maxiter)

logger.info(f"total it time={avg_time}")

avg_time /= iters

logger.info(f"avg it time={avg_time}")

if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

error = HDG_wrapper.errors(vectorSolution)[0]
logger.info(f"HDG_wrapper error={error:>.6e}")

HDG_wrapper.plot_option("outputDir", args.output_dir)
HDG_wrapper.plot_option("fileName", args.output)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_solution(vectorSolution)

logger.info(f"solution written to '{output_path}'")
