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

import prin2, jprecond

# '<' for little-endian
# u8 : uint64_t
# S8 : char[8]

DataTableType = np.dtype([
  ('name', 'S8'),
  ('offset', '<u8'),
  ('size', '<u8'),
])

GeoBinHeaderType = np.dtype([
    ('magic', 'S8'),
    ('space_dim', '<u8'),
    ('hyperedge_dim', '<u8'),
    ('n_points', '<u8'),
    ('n_hypernodes', '<u8'),
    ('n_hyperedges', '<u8'),
])

FloatType = np.dtype("float64")

######### SETUP argument parsing & logging

parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("network")
parser.add_argument("-t", "--rtol",
  help="relative tolerance when to stop the CG iterator",
  type=float, default=1e-10
)
parser.add_argument("--txt", help="expect .txt instead of .bin files", action="store_true")
parser.add_argument("-d", "--debug", help="toggle debug mode", action="store_true")

logging.setLoggerClass(prin2.Logger)
logger = logging.getLogger("fiber_network_elastic")

args = parser.parse_args()
logger.log_args(args)
os.system("mkdir -p output")

network_path = args.network + ".geo"
if not args.txt:
  network_path += ".bin"
logger.info(f"{network_path=}")

######## MAIN code

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG
  
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.local_solver    = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = args.debug

logger.info("reading files")

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP(network_path)
if args.txt:
  network_points = np.loadtxt(args.network + ".pts")
else:
  header_without_tables = np.fromfile(
    network_path,
    dtype=GeoBinHeaderType,
    offset=0,
    count=1
  )[0]
  tables = np.fromfile(
    network_path,
    dtype=DataTableType,
    offset=GeoBinHeaderType.itemsize,
    count=5
  )
  network_points = np.fromfile(
    network_path,
    dtype=FloatType,
    offset=tables[0]["offset"],
    count=int(tables[0]["size"]/FloatType.itemsize),
  ).reshape(-1, header_without_tables["space_dim"])

logger.info("computing residual")

rhs = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

logger.info("assembling  A...")

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

logger.info("assembling  B...")

B = sp.linalg.LinearOperator(
  (system_size,system_size),
  matvec=jprecond.JPrecond(A, network_points, [2**3, 2**3], repeat=6).matmul
)

iters = 0
start = datetime.datetime.now()

def log_iter(x):
  global iters, start
  iters += 1
  relErr = np.linalg.norm(A @ x - rhs) / np.linalg.norm(rhs)
  bilin = .5 * x.dot(A @ x) - x.dot(rhs)
  duration = datetime.datetime.now() - start
  start = datetime.datetime.now()
  logger.info(f"{iters:>5} {relErr:>13.6e} {bilin:>13.6e} {duration}")

logger.info("starting cg")

vectorSolution, num_iter = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=log_iter, M=B)

if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

error = HDG_wrapper.errors(vectorSolution)[0]
logger.info(f"HDG_wrapper error={error:>.6e}")

output_name = args.network + "_timo"
HDG_wrapper.plot_option("fileName", output_name)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution)

logger.info(f"solution written to 'output/{output_name}'")
