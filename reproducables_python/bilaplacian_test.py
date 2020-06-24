# Python running example for HDG solution of an elasticity problem on a superaggregate!

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
problem = "AbstractProblem \
  < \
    Topology::Cubic<1,1>, Geometry::UnitCube<1,1,double>, \
    NodeDescriptor::Cubic<1,1>, \
    bilaplacian<1,1,2,TestParametersSin,double> \
  >"
filenames = [ "HyperHDG/Geometry/Cubic.hxx" , \
              "HyperHDG/NodeDescriptor/Cubic.hxx", \
              "HyperHDG/LocalSolver/BernoulliBeams.hxx", \
              "reproducables_python/parameters/bilaplacian.hxx" ]

# Import C++ wrapper class to use HDG method on graphs.
from cython_import import cython_import
PyDP = cython_import(["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames, True)

# Initialising the wrapped C++ class HDG_wrapper.
HDG_wrapper = PyDP( [4] )

# Generate right-hand side vector "vectorRHS = - A * vectorDirichlet", where vectorDirichlet is the
# vector of Dirichlet values.
vectorRHS = HDG_wrapper.return_zero_vector();
vectorRHS = [-i for i in HDG_wrapper.total_flux_vector(vectorRHS)]

# Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
# fashion.
system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

# Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
[vectorSolution, num_iter] = sp_lin_alg.bicgstab(A, vectorRHS, tol=1e-9) # Parameters for BiCGStab.

# Check, whether solution has been successful!
if num_iter == 0:
  print("Solution has been successfully calculated!")
else:
  print("The linear solver (conjugate gradients) failed (did not converge)!")

print("Error: " + str(HDG_wrapper.calculate_L2_error(vectorSolution)))
