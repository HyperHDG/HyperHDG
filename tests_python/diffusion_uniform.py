from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

import os, sys

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG

const                     = HyperHDG.config()
const.global_loop         = "Elliptic"
const.local_solver        = "DiffusionUniform < 1, 1, 2 * 1 >"
const.topology            = "Cubic< 1, 3 >"
const.geometry            = "UnitCube< 1, 3 >"
const.node_descriptor     = "Cubic< 1, 3 >"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.debug_mode          = True
const.allow_file_output   = False

PyDiffusionProblem = HyperHDG.include(const)

tolerance   = 1e-8
HDG_wrapper = PyDiffusionProblem([4,2,2])

vectorDirichlet    = HDG_wrapper.zero_vector()
vectorDirichlet[0] = 1.

index_vector = np.array([ 0, len(vectorDirichlet)-1 ])
HDG_wrapper.read_dirichlet_indices(index_vector)

vectorRHS = [-i for i in HDG_wrapper.trace_to_flux(vectorDirichlet)]

system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.trace_to_flux )

[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=100, tol=1e-9) # Parameters for CG.

reference_solution = np.array(
  [ 1.,         0.6999695,  0.55280737, 0.46359316, 0.41591649, 0.72849089,
    0.62353531, 0.52383342, 0.4428244,  0.39207816, 0.63876017, 0.57986777,
    0.5,        0.42013223, 0.36123983, 0.72849089, 0.62353531, 0.52383342,
    0.4428244,  0.39207816, 0.65166809, 0.58551499, 0.5,        0.41448501,
    0.34833191, 0.60792184, 0.5571756,  0.47616658, 0.37646469, 0.27150911,
    0.63876017, 0.57986777, 0.5,        0.42013223, 0.36123983, 0.60792184,
    0.5571756,  0.47616658, 0.37646469, 0.27150911, 0.58408351, 0.53640684,
    0.44719263, 0.3000305,  0. ])

assert np.linalg.norm(reference_solution - vectorSolution - vectorDirichlet, np.inf) < tolerance
