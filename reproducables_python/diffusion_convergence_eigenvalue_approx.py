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
    self.hdg_wrapper       = hdg_wrapper
    self.dirichlet_inidces = hdg_wrapper.dirichlet_indices()
  def long_vector_size(self):
    return self.hdg_wrapper.size_of_system()
  def short_vector_size(self):
    return self.hdg_wrapper.size_of_system() - len(self.dirichlet_inidces)
  def long_vector(self, vector):
    assert len(vector) == self.short_vector_size()
    vec = [0.] * (len(vector)+len(self.dirichlet_inidces))
    n_indices = 0
    for i in range(len(vector)):
      while n_indices < len(self.dirichlet_inidces) and \
            i + n_indices == self.dirichlet_inidces[n_indices]:
        n_indices = n_indices + 1
      vec[i + n_indices] = vector[i]
    return vec
  def short_vector(self, vector):
    assert len(vector) == self.long_vector_size()
    vec = [0.] * (len(vector)-len(self.dirichlet_inidces))
    n_indices = 0
    for i in range(len(vec)):
      while n_indices < len(self.dirichlet_inidces) and \
            i + n_indices == self.dirichlet_inidces[n_indices]:
        n_indices = n_indices + 1
      vec[i] = vector[i + n_indices]
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
  def multiply_inverse_mass(self, vector):
    system_size = self.short_vector_size()
    Mass        = LinearOperator( (system_size,system_size), matvec= self.multiply_mass )
    [vec, n_it] = sp_lin_alg.cg(Mass, vector, tol=1e-9)
    assert n_it == 0
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx(poly_degree, dimension, iteration):
  
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Diffusion<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersEigs,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/Diffusion.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )
  helper = helper_ev_approx(HDG_wrapper)
  
  # Define LinearOperator in terms of C++ functions to use scipy's matrix-free linear solvers.
  system_size = helper.short_vector_size()
  Stiff = LinearOperator( (system_size,system_size), matvec= helper.multiply_stiff )
  Mass  = LinearOperator( (system_size,system_size), matvec= helper.multiply_mass )
  M_inv = LinearOperator( (system_size,system_size), matvec= helper.multiply_inverse_mass )
  
  # vector = [0.] * system_size
  # for i in range(system_size):
  #   vector[i] = 1.
  #   print(Mass * vector)
  #   vector[i] = 0.
  
  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
  [vals, vecs] = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, which='SM', tol=1e-9, Minv=M_inv)

  # Print error.
  error = np.absolute(vals[0] - np.pi * np.pi)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_convergence_eigenvalue_approx.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "diff_e-" + str(dimension) + "-" + str(iteration) );
  HDG_wrapper.plot_option( "printFileNumber" , "false" );
  HDG_wrapper.plot_option( "scale" , "0.95" );
  HDG_wrapper.plot_solution(helper.long_vector([x[0].real for x in vecs]));
  
  # Return smallest eigenvalue and corresponding eigenvector.
  return vals[0].real, helper.long_vector([x[0].real for x in vecs])
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(2,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(2,3):
        eigenvalue_approx(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
