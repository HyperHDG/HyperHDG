# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

# Import package to print date and time of start and end of program.
from datetime import datetime

# Correct the python paths!
import os, sys
sys.path.append(os.path.dirname(__file__) + "/..")

# --------------------------------------------------------------------------------------------------
# THIS SECTION CAN BE CHANGED:

# Define aggregate specification:
aggregate = "5"
# aggregate = "1000_tree"
# aggregate = "5000_tree"
# --------------------------------------------------------------------------------------------------

# Print starting time of diffusion test.
start_time = datetime.now()
print("Starting time is", start_time)

# Predefine problem to be solved.
problem = "GlobalLoop::Elliptic< Topology::File<1,3>, Geometry::File<1,3>, NodeDescriptor::File<1,3>, "\
         +                  "LocalSolver::LengtheningBernoulliBendingBeam<1,3,1,2> > "
filenames = [ "HyperHDG/geometry/file.hxx" , \
              "HyperHDG/node_descriptor/file.hxx", \
              "HyperHDG/local_solver/bernoulli_beams.hxx" ]

# Import C++ wrapper class to use HDG method on graphs.
from cython_import import cython_import
PyDP = cython_import(["elliptic_loop", problem, "string", "string"], filenames, True)

# Initialising the wrapped C++ class HDG_wrapper.
HDG_wrapper = PyDP( "domains/aggregate_" + aggregate + ".pts" )

# Initialize vector containing the Dirichlet values: Indices not set in the index_vector are ignored
# here. However, values not equal zero in vectorDirichlet that have indices that do not occur in the
# index vector (next) will cause a wrong representation of the final result.
vectorDirichlet = HDG_wrapper.return_zero_vector()
vectorDirichlet[0] = 1.

# Set the hypernodes that are supposed to be of Dirichlet type.
# Note that all non-zero entries of vectorDirichlet are supposed to be contained in the index vector
# to keep consistency.
index_vector = np.array([ 0, 1, 2, 3, \
  len(vectorDirichlet)-4, len(vectorDirichlet)-3, len(vectorDirichlet)-2, len(vectorDirichlet)-1 ])
HDG_wrapper.read_dirichlet_indices(index_vector)

# Generate right-hand side vector "vectorRHS = - A * vectorDirichlet", where vectorDirichlet is the
# vector of Dirichlet values.
vectorRHS = [-i for i in HDG_wrapper.matrix_vector_multiply(vectorDirichlet)]

# Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
# fashion.
system_size = HDG_wrapper.size_of_system()
A = LinearOperator( (system_size,system_size), matvec= HDG_wrapper.matrix_vector_multiply )

# Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, maxiter=5500, tol=1e-9) # Parameters for CG.

# Check, whether solution has been successful!
if num_iter == 0:
  print("Solution has been successfully calculated!")
else:
  print("The linear solver (conjugate gradients) failed (did not converge)!")

# Plot solution to vtu File to be visualized using Paraview.
HDG_wrapper.plot_option("fileName", "aggregate_" + aggregate)
HDG_wrapper.plot_option("printFileNumber", "false" )
HDG_wrapper.plot_option("plotEdgeBoundaries", "true")
HDG_wrapper.plot_option("scale", "0.8")
HDG_wrapper.plot_option("boundaryScale", "0.9")
HDG_wrapper.plot_solution(vectorSolution + vectorDirichlet)
print("Solution written to file" , HDG_wrapper.plot_option("fileName", ""), "in output directory.")

# Print ending time of diffusion test.
end_time = datetime.now()
print("Program ended at", end_time, "after", end_time-start_time)