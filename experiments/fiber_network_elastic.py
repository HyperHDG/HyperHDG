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

DomainsHeaderType = np.dtype([
  ('magic', 'S8'),
  ('idsize', '<u8'),
  ('n_domains', '<u8'),
  ('tables', DataTableType, (2,))
])

FloatType = np.dtype("float64")
IdType = np.dtype("<u4")

######### SETUP argument parsing & logging

default_output_dir = "output"
now = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
default_output_name = f"{os.path.basename(__file__)}.{now}"
parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("network", help="full path to the network file")
parser.add_argument("-t", "--rtol",
  help="relative tolerance when to stop the CG iterator",
  type=float, default=1e-10
)
parser.add_argument("-d", "--debug", help="toggle debug mode", action="store_true")
parser.add_argument("-o", "--output", help="output name", default=default_output_name)
parser.add_argument("--output-dir", help="output dir", default=default_output_dir)
parser.add_argument("--domains", help="domains file")
parser.add_argument("--maxiter", help="maximum number of cg iterations", type=int, default=100)

logging.setLoggerClass(prin2.Logger)
logger = logging.getLogger("fiber_network_elastic")

args = parser.parse_args()
logger.log_args(args)

# verify output path is writable
output_path = f"{args.output_dir}/{args.output}"
try:
  with open(output_path, "w") as file:
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
const.local_solver    = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = args.debug

logger.info("reading files")

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP(args.network)

with zstd.open(args.network, "rb") as decom_file:
  header_without_tables = np.frombuffer(
    decom_file.read(GeoBinHeaderType.itemsize),
    dtype=GeoBinHeaderType,
    offset=0,
    count=1
  )[0]
  tables = np.frombuffer(
    decom_file.read(5*DataTableType.itemsize),
    dtype=DataTableType,
    count=5
  )
  size = tables[0]["size"]
  network_points = np.frombuffer(
    decom_file.read(size),
    dtype=FloatType,
    count=int(size/FloatType.itemsize),
  ).reshape(-1, header_without_tables["space_dim"])

domains = None
if args.domains:
  with zstd.open(args.domains, "rb") as file:
    header = np.frombuffer(
      file.read(DomainsHeaderType.itemsize),
      dtype=DomainsHeaderType,
    )[0]
    ioffsets = np.frombuffer(
      file.read(header["tables"][0]["size"]),
      dtype=IdType,
    )
    all_domains = np.frombuffer(
      file.read(header["tables"][1]["size"]),
      dtype=IdType,
    )
    domains = jprecond.Domains(ioffsets, all_domains)

logger.info("computing residual")

rhs = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

logger.info("assembling  A...")

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csc_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
logger.info(f"{A.shape=}")

logger.info("assembling  B...")

B = sp.linalg.LinearOperator(
  (system_size,system_size),
  # repeat=6 because for every node we have 6 unknowns, namely displacement+rotation
  matvec=jprecond.JPrecond(A, network_points, [2**3, 2**3], repeat=6, domains=domains).matmul
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

vectorSolution, num_iter = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=log_iter, M=B, maxiter=args.maxiter)

if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

error = HDG_wrapper.errors(vectorSolution)[0]
logger.info(f"HDG_wrapper error={error:>.6e}")

HDG_wrapper.plot_option("outputDir", args.output_dir)
HDG_wrapper.plot_option("fileName", args.output)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution)

logger.info(f"solution written to '{output_path}'")
