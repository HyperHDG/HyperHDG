from __future__ import print_function

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

import matplotlib.pyplot as plt

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
const.local_solver    = "LengtheningBeam<1,2,1,2>"
# const.local_solver    = "Bilaplacian_cont_flux<1,4,8>"
# const.local_solver    = "Bilaplacian<1,1,2>"
# const.topology        = "File<1,3>"
# const.geometry        = "File<1,3>"
# const.node_descriptor = "File<1,3>"
const.topology        = "Cubic<1,1>"
const.geometry        = "UnitCube<1,1,double>"
const.node_descriptor = "Cubic<1,1>"
# const.cython_replacements = ["string", "string"]
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.include_files   = [ "HyperHDG/local_solver/diffusion_ldgh.hxx", \
                          "HyperHDG/local_solver/bilaplacian_ldgh_cont_flux.hxx" ]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
# HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
#                     "/../domains/aggregate_" + aggregate + ".pts" )
# HDG_wrapper = PyDP( os.path.dirname(os.path.abspath(__file__)) + \
#     "/../domains/fiber_network_1000.geo" )
HDG_wrapper = PyDP( [4] * 1 )

vectorDirichlet = HDG_wrapper.zero_vector()
# vectorDirichlet[0] = 1.

# index_vector = np.array([ 0, 4 ])
# HDG_wrapper.read_dirichlet_indices(index_vector)

vectorRHS = [-i for i in HDG_wrapper.residual_flux(vectorDirichlet)]
# print(vectorRHS)

system_size = HDG_wrapper.size_of_system()
# A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

# B = A.todense()
# B = B[2:-2,2:-2]
# print(B)

# col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat()
# C = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))

# print(np.subtract(A.dot(np.ones(system_size)),C.dot(np.ones(system_size))))

vectorSolution, num_iter = sp_lin_alg.cg(A, vectorRHS, tol=1e-9)
if num_iter != 0:
  raise RuntimeError("Linear solver did not converge!")

# print("The linear solver (conjugate gradients) worked!")



print(vectorSolution)

# vectorSolution = np.array([0,0,.5,1,0,0])
# print(vectorSolution)
print(A.dot(vectorSolution))
print(vectorRHS)
print(HDG_wrapper.residual_flux(vectorSolution))



# print(vectorRHS)

error = HDG_wrapper.errors(vectorSolution)[0]
print(" Error: ", error)

# print(vectorSolution)

HDG_wrapper.plot_option("fileName", "aggregate_" + aggregate)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution + vectorDirichlet)
print("Solution written to file" , HDG_wrapper.plot_option("fileName", ""), "in output directory.")

end_time = datetime.now()
print("Program ended at", end_time, "after", end_time-start_time)
