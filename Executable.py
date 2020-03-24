# Python running example for HDG solution of a diffusion problem which is formulated on a graph!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Predefine problem to be solved.
problem = "AbstractProblem < Topology::Cubic< 1, 3 >, " \
         +                  "Geometry::UnitCube< 1, 3 >, " \
         +                  "Diffusion_TensorialUniform < 1, 1, 2 * 1 > " \
         +                ">"

# Import C++ wrapper class to use HDG method on graphs.
from hyImport import hyImport
PyDP = hyImport(["AbstractProblem", problem])

# Initialising the wrapped C++ class HDG_wrapper.
HDG_wrapper = PyDP([1,1,1])
# HDG_wrapper = PyElasticityProblem("domains/SimpleTriangle.geo")

# Initialize vector containing the Dirichlet values: Indices not set in the index_vector are ignored
# here. However, values not equal zero in vectorDirichlet that have indices that do not occur in the
# index vector (next) will cause a wrong representation of the final result.
vectorDirichlet = HDG_wrapper.return_zero_vector()
vectorDirichlet[0] = 1.
vectorDirichlet[2] = 2.
vectorDirichlet[4] = 3.
vectorDirichlet[6] = 4.
# vectorDirichlet[len(vectorDirichlet)-1] = 1. # Comment if checking for trivial solution.

# Set the hypernodes that are supposed to be of Dirichlet type.
# Note that all non-zero entries of vectorDirichlet are supposed to be contained in the index vector
# to keep consistency.
index_vector = np.array([ 0, 1, 2, 3, 4, 5, 6, 7])# len(vectorDirichlet)-1 ])
HDG_wrapper.read_dirichlet_indices(index_vector)

# Print index vector and vector containing the Dirichlet values.
print("Dirichlet indices: ", index_vector)
print("Dirichlet values: ", vectorDirichlet)

# Generate right-hand side vector "vectorRHS = - A * vectorDirichlet", where vectorDirichlet is the
# vector of Dirichlet values.
vectorRHS = [-i for i in HDG_wrapper.matrix_vector_multiply(vectorDirichlet)]

# Print right-hand side vector.
print("Right-hand side: ", vectorRHS)

# Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
# fashion.
system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

# Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=100, tol=1e-9) # Parameters for CG.

# Print Solution to the problem (which is x + x_D, i.e. vectorSolution + vectorDirichlet) or number
# of CG iterations num_iter.
if num_iter == 0:
  print("Solution:\n", vectorSolution + vectorDirichlet)
else:
  print("The linear solver (conjugate gradients) failed (did not converge)!")

# Plot solution to vtu File to be visualized using Paraview.
HDG_wrapper.plot_solution(vectorSolution + vectorDirichlet)
print("Solution written to file" , HDG_wrapper.plot_option("fileName", ""), "in output directory.")
