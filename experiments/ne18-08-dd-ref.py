# --------------------------------------------------------------------------------------------------
# The data necessary to run this file can be obtained using the HyperHDG.fiber_network.make_geo
# function on the files of Hauck, M., & Rupp, A. (2024). Fiber network models of paper. Zenodo.
# https://doi.org/10.5281/zenodo.12751486
# --------------------------------------------------------------------------------------------------

from __future__ import print_function

import numpy as np
import scipy.sparse as sp

from datetime import datetime
import time

import os, sys
import argparse

# --------------------------------------------------------------------------------------------------
# THIS SECTION CAN BE CHANGED:

# --------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
  prog='ne18-08-dd-ref',
  description='solve fiber network with the reference implementation',
)
parser.add_argument("--domain")
parser.add_argument("--points")
parser.add_argument("--subdomains", type=int)
args = parser.parse_args()
domain = args.domain

start_time = datetime.now()
print("Starting time is", start_time)
os.system("mkdir -p output")

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG

const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.local_solver    = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoshenkoClampedConstant>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP(domain)

vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )


print("Start matrix setup:", datetime.now())

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

print("A setup", datetime.now())

points = np.loadtxt(args.points)
helper = HyperHDG.fiber_network.precond(points, [args.subdomains, args.subdomains], repeat=6)
def precond_mult( vec_x ):
  return helper.precond(A, vec_x)
B = sp.linalg.LinearOperator( (system_size,system_size), matvec= precond_mult )

print("B setup", datetime.now())

# Emit the same YAML block as the C++ solvers' -ksp_monitor_yaml (KSPMonitorYAML
# in prin2.cxx): one "- it/time/rnorm" entry per iteration, rnorm the absolute
# unpreconditioned residual ||b - A x||_2 (matching -ksp_norm_type unpreconditioned),
# time in seconds since the start of the solve.
iters = 0
ksp_t0 = 0.
def nonlocal_iterate(vec_x):
  global iters
  iters += 1
  print("  - it: %3d" % iters)
  print("    time: %.16e" % (time.perf_counter() - ksp_t0))
  print("    rnorm: %.16e" % np.linalg.norm(vectorRHS - A.dot(vec_x)))

print("ksp_monitor:")
ksp_t0 = time.perf_counter()
# it 0: zero initial guess -> r0 = b, so rnorm = ||b||
print("  - it: %3d" % 0)
print("    time: %.16e" % (time.perf_counter() - ksp_t0))
print("    rnorm: %.16e" % np.linalg.norm(vectorRHS))

vectorSolution, num_iter = sp.linalg.cg(A, vectorRHS, rtol=1e-10, callback=nonlocal_iterate, M=B)
if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")


error = HDG_wrapper.errors(vectorSolution)[0]
print(" Error: ", error)

# print(vectorSolution)

HDG_wrapper.plot_option("fileName", domain + "_timo")
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution)
print("Solution written to file" , HDG_wrapper.plot_option("fileName", ""), "in output directory.")

end_time = datetime.now()
print("Program ended at", end_time, "after", end_time-start_time)
