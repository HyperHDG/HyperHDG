# --------------------------------------------------------------------------------------------------
# The data necessary to run this file can be obtained using the HyperHDG.fiber_network.make_geo
# function on the files of Hauck, M., & Rupp, A. (2024). Fiber network models of paper. Zenodo.
# https://doi.org/10.5281/zenodo.12751486
# --------------------------------------------------------------------------------------------------

from __future__ import print_function

import numpy as np
import scipy.sparse as sp

from datetime import datetime

import os, sys

# --------------------------------------------------------------------------------------------------
# THIS SECTION CAN BE CHANGED:

# domain = "fiber_network_1000"
# domain = "fiber_network_14871"
domain = "fiber_network_615452"
# --------------------------------------------------------------------------------------------------

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
const.local_solver    = "TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + "/../domains/" + domain + ".geo" )

vectorRHS = np.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )


print("Start matrix setup:", datetime.now())

system_size = HDG_wrapper.size_of_system()
col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

print("A setup", datetime.now())

points = np.loadtxt("domains/" + domain + "_points.txt")
helper = HyperHDG.fiber_network.precond(points, [2**3, 2**3], repeat=6)
def precond_mult( vec_x ):
  return helper.precond(A, vec_x)
B = sp.linalg.LinearOperator( (system_size,system_size), matvec= precond_mult )

print("B setup", datetime.now())

iters = 0
def nonlocal_iterate(vec_x):
  global iters
  iters += 1
  print(iters, "\t", np.linalg.norm(A.dot(vec_x) - vectorRHS) / np.linalg.norm(vectorRHS),
        "\t", .5 * vec_x.dot(A.dot(vec_x)) - vec_x.dot(vectorRHS), " \t", datetime.now())

print(np.linalg.norm(vectorRHS))

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
