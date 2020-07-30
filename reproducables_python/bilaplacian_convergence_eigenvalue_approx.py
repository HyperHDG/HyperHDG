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


# --------------------------------------------------------------------------------------------------
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_approx():
  def __init__(self, hdg_wrapper):
    self.hdg_wrapper = hdg_wrapper
  def long_vector(self, vector):
    vec = [0.] * (len(vector)+4)
    for i in range(len(vector)):
      vec[i+2] = vector[i]
    return vec
  def short_vector(self, vector):
    vec = [0.] * (len(vector)-4)
    for i in range(len(vec)):
      vec[i] = vector[i+2]
    return vec
  def multiply_stiff(self, vector):
    vec = self.long_vector(vector)
    vec = np.multiply(self.hdg_wrapper.matrix_vector_multiply(vec), -1.)
    vec = self.short_vector(vec)
    return vec
  def multiply_mass(self, vector):
    vec = self.long_vector(vector)
    vec = self.hdg_wrapper.mass_matrix_multiply(vec)
    vec = self.short_vector(vec)
    return vec
  def plot_vector(self, vector):
    vec = self.long_vector(vector)
    for i in range(len(vec)):
      vec[i] = (vec[i]).real
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx(poly_degree, dimension, iteration):
  
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "bilaplacian<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersEigs,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/bilaplacian.hxx", \
                "reproducables_python/parameters/bilaplacian.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )
  helper = helper_ev_approx(HDG_wrapper)

  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system()
  Stiff = LinearOperator( (system_size-4,system_size-4), matvec= helper.multiply_stiff )
  Mass  = LinearOperator( (system_size-4,system_size-4), matvec= helper.multiply_mass )
  
  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
  [vals, vecs] = sp_lin_alg.eigs(Stiff, k=1, M=Mass, which='SM', tol=1e-9)

  # Print error.
  error = np.absolute(vals[0] - np.power(np.pi, 4))
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/bilaplacian_convergence_eigenvalue_approx.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "diff_e-" + str(dimension) + "-" + str(iteration) );
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(helper.plot_vector(vecs));
  
  return vals[0], helper.long_vector(vecs)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,2):
      print("Dimension is ", dimension, "\n")
      for iteration in range(2,8):
        eigenvalue_approx(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
