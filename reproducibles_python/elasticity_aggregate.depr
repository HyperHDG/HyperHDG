from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

from datetime import datetime

import os, sys

# --------------------------------------------------------------------------------------------------
# THIS SECTION CAN BE CHANGED:

# Define aggregate specification:
aggregate = "5"
# aggregate = "1000_tree"
# aggregate = "5000_tree"
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
const.local_solver    = "LengtheningBernoulliBendingBeam<1,3,1,2>"
const.topology        = "File<1,3>"
const.geometry        = "File<1,3>"
const.node_descriptor = "File<1,3>"
const.cython_replacements = ["string", "string"]
const.include_files   = [ "HyperHDG/local_solver/diffusion_uniform_ldgh.hxx", \
                          "HyperHDG/local_solver/bilaplacian_uniform_ldgh.hxx" ]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
                    "/../domains/aggregate_" + aggregate + ".pts" )

vectorDirichlet = HDG_wrapper.zero_vector()
vectorDirichlet[0] = 1.

index_vector = np.array([ 0, 1, 2, 3, \
  len(vectorDirichlet)-4, len(vectorDirichlet)-3, len(vectorDirichlet)-2, len(vectorDirichlet)-1 ])
HDG_wrapper.read_dirichlet_indices(index_vector)

vectorRHS = [-i for i in HDG_wrapper.trace_to_flux(vectorDirichlet)]

system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=5500, tol=1e-9)
if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

print("The linear solver (conjugate gradients) worked!")

HDG_wrapper.plot_option("fileName", "aggregate_" + aggregate)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution + vectorDirichlet)
print("Solution written to file" , HDG_wrapper.plot_option("fileName", ""), "in output directory.")

end_time = datetime.now()
print("Program ended at", end_time, "after", end_time-start_time)
