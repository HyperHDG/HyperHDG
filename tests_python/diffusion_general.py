# Python running example for HDG solution of a diffusion problem which is formulated on a graph!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Correct the python paths!
import os, sys
sys.path.append(os.path.dirname(__file__) + "/..")


# Predefine problem to be solved.
problem = "AbstractProblem < Topology::Cubic< 1, 3 >, " \
         +                  "Geometry::UnitCube< 1, 3 >, " \
         +                  "NodeDescriptor::Cubic< 1, 3 >, " \
         +                  "Diffusion< 1, 3, 1, 2 * 1 > " \
         +                ">"
filenames = [ "HyperHDG/Geometry/Cubic.hxx" , \
              "HyperHDG/NodeDescriptor/Cubic.hxx" , \
              "HyperHDG/LocalSolver/Diffusion.hxx" ]

# Import C++ wrapper class to use HDG method on graphs.
from cython_import import cython_import
PyDiffusionProblem = cython_import \
  (["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames)

# Define tolerance
tolerance = 1e-8

# Initialising the wrapped C++ class HDG_wrapper.
HDG_wrapper = PyDiffusionProblem([4,2,2])

# Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
# fashion.
system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

# Print index vector and vector containing the Dirichlet values.
# print("Dirichlet indices: ", index_vector)
# print("Dirichlet values: ", vectorDirichlet)

# Generate right-hand side vector "vectorRHS = - A * vectorDirichlet", where vectorDirichlet is the
# vector of Dirichlet values.
vectorRHS = HDG_wrapper.return_zero_vector()
vectorRHS = [-i for i in HDG_wrapper.total_flux_vector(vectorRHS)]

# Print right-hand side vector.
print("Right-hand side:\n", vectorRHS)



# Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=100, tol=1e-9) # Parameters for CG.

# Print Solution to the problem (which is x + x_D, i.e. vectorSolution + vectorDirichlet) or number
# of CG iterations num_iter.
if num_iter == 0:
  print("Solution:\n", vectorSolution)
else:
  print("The linear solver (conjugate gradients) failed after", num_iter, " iterations.")

# Plot solution to vtu File to be visualized using Paraview.
HDG_wrapper.plot_solution(vectorSolution)
print("Solution written to file", HDG_wrapper.plot_option("fileName", ""), "in output directory.")

reference_solution = np.array(
  [ 0.,         0.6999695,  0.55280737, 0.46359316, 0.41591649, 0.72849089,
    0.62353531, 0.52383342, 0.4428244,  0.39207816, 0.63876017, 0.57986777,
    0.5,        0.42013223, 0.36123983, 0.72849089, 0.62353531, 0.52383342,
    0.4428244,  0.39207816, 0.65166809, 0.58551499, 0.5,        0.41448501,
    0.34833191, 0.60792184, 0.5571756,  0.47616658, 0.37646469, 0.27150911,
    0.63876017, 0.57986777, 0.5,        0.42013223, 0.36123983, 0.60792184,
    0.5571756,  0.47616658, 0.37646469, 0.27150911, 0.58408351, 0.53640684,
    0.44719263, 0.3000305,  0. ])

# if np.linalg.norm(reference_solution - vectorSolution, np.inf) > tolerance:
#   print("Diffusion test FAILED!")
# else:
#   print("Diffusion test SUCCEEDED!")
