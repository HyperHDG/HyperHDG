# --------------------------------------------------------------------------------------------------
# The data necessary to run this file can be obtained using the HyperHDG.fiber_network.make_geo
# function on the files of Hauck, M., & Rupp, A. (2024). Fiber network models of paper. Zenodo.
# https://doi.org/10.5281/zenodo.12751486
# --------------------------------------------------------------------------------------------------

from __future__ import print_function
import numpy as np
import scipy.sparse as sp
import os, sys, argparse, logging, datetime, subprocess

# SETUP argument parsing & logging

time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
git_hash =  subprocess.run(["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()

logger = logging.getLogger("fiber_network_elastic")
logger.setLevel(logging.INFO)
log_formatter  = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

fhandler = logging.FileHandler(f"fiber_network_elastic.{time_stamp}.log", mode="w")
fhandler.setFormatter(log_formatter)
fhandler.setLevel(logging.INFO)

shandler = logging.StreamHandler(sys.stdout)
shandler.setFormatter(log_formatter)
shandler.setLevel(logging.INFO)

logger.addHandler(fhandler)
logger.addHandler(shandler)

parser = argparse.ArgumentParser(description="fiber_network_elastic by Joseph Holten")
parser.add_argument("network")
parser.add_argument("-t", "--rtol",
  help="relative tolerance when to stop the CG iterator",
  type=float, default=1e-10
)

args = parser.parse_args()
domain = args.network

# MAIN code

logger.info("started with")
logger.info(f"  git_hash={git_hash}")
logger.info(f"  network={args.network}")
logger.info(f"  rtol={args.rtol}")

os.system("mkdir -p output")

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
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP( domain + ".geo" )

rhs = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

logger.info("assembling  A...")

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

logger.info("assembling  B...")

points = np.loadtxt(domain + "_points.txt")
helper = HyperHDG.fiber_network.precond(points, [2**3, 2**3], repeat=6)

B = sp.linalg.LinearOperator( (system_size,system_size), matvec= lambda x: helper.precond(A, x))

iters = 0

def log_iter(x):
  global iters
  iters += 1
  relErr = np.linalg.norm(A @ x - rhs) / np.linalg.norm(rhs)
  bilin = .5 * x.dot(A @ x) - x.dot(rhs)
  logger.info(f"{iters:>5} {relErr:>13.6e} {bilin:>13.6e}")

logger.info("starting cg")

vectorSolution, num_iter = sp.linalg.cg(A, rhs, rtol=args.rtol, callback=log_iter, M=B)

if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

error = HDG_wrapper.errors(vectorSolution)[0]
logger.info(f"HDG_wrapper error={error:>.6e}")

output_name = domain + "_timo"
HDG_wrapper.plot_option("fileName", output_name)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution)

logger.info(f"solution written to 'output/{output_name}'")
