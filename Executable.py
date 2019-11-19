# Python running example for HDG solution of a diffusion problem which is
# formulated on a graph!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths
# libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Import C++ wrapper class to use HDG method on graphs.
from ClassWrapper import PyDiffusionProblem


# Initialising the wrapped c++ class HDG_wrapper.
HDG_wrapper = PyDiffusionProblem([2,2])

# Set the hypernodes that are supposed to be of Dirichlet type.
index_vector = np.array([ 0, 8 ])
HDG_wrapper.read_dirichlet_indices(index_vector)

# Initialize vector containing the Dirichlet values: Indices not set in the
# index_vector are ignored here. However, values != 0 in vectorDirichlet
# that have indices that do not occur in the index vector will cause a wrong
# representation of the final result.
vectorDirichlet = HDG_wrapper.return_zero_vector()
vectorDirichlet[0] = 1.
vectorDirichlet[8] = 0.

# Print index vector and vector containing the Dirichlet values.
print("Dirichlet indices: ", index_vector)
print("Dirichlet values: ", vectorDirichlet)

# Generate right-hand side vector "b = - A * x_D", where x_D is the vector
# of Dirichlet values.
vectorRHS = HDG_wrapper.matrix_vector_multiply(vectorDirichlet)
vectorRHS = [-i for i in vectorRHS]

# Print right-hand side vector.
print("Right-hand side: ", vectorRHS)

# Define LinearOperator to use scipy linear solvers in a matrix-free fashion.
system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size),
                    matvec= HDG_wrapper.matrix_vector_multiply )

# Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=100)

# Print Solution to the problem (which is x + x_D) or number of CG iterations.
if num_iter == 0:
  print("Solution: ", vectorSolution + vectorDirichlet)
else:
  print("CG failde with total number of ", num_iter, " iterations.")
