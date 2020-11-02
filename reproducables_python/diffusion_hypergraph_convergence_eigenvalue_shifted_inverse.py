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
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_approx():
  def __init__(self, hdg_wrapper, sigma):
    self.hdg_wrapper       = hdg_wrapper
    self.sigma             = sigma
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
    vec = self.hdg_wrapper.matrix_vector_multiply(self.long_vector(vector))
    vec = np.multiply(self.short_vector(vec), -1.)
    return vec
  def shifted_mult(self, vector):
    vec = self.hdg_wrapper.matrix_vector_multiply(self.long_vector(vector), self.sigma)
    vec = np.multiply(self.short_vector(vec), -1.)
    return vec
  def multiply_mass(self, vector):
    vec = self.hdg_wrapper.matrix_vector_multiply(self.long_vector(vector), self.sigma)
    vec = np.subtract( self.hdg_wrapper.matrix_vector_multiply(self.long_vector(vector)) , vec )
    vec = np.multiply( self.short_vector(vec) , -1. / self.sigma )
    return vec
  def shifted_inverse(self, vector):
    if not hasattr(self, 'ShiftOp'):
      system_size  = self.short_vector_size()
      self.ShiftOp = LinearOperator( (system_size,system_size), matvec=self.shifted_mult )
    [vec, n_it] = sp_lin_alg.cg(self.ShiftOp, vector, tol=1e-9)
    assert n_it == 0
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def eigenvalue_approx_SI(poly_degree, dimension, iteration):
  # Print starting time of diffusion test.
  start_time = datetime.now()
  print("Starting time is", start_time)
  
  # Predefine problem to be solved.
  problem = "GlobalLoop::ShiftedEigenvalue < Topology::Cubic<" + str(dimension) + ",3>, " \
          + "Geometry::UnitCube<" + str(dimension) + ",3,double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + ",3>, " \
          + "LocalSolver::DiffusionEigs<" + str(dimension) + "," + str(poly_degree) + "," \
          + str(2*poly_degree) + ",TestParametersEigs,double> >"
  filenames = [ "HyperHDG/geometry/cubic.hxx" , "HyperHDG/node_descriptor/cubic.hxx", \
                "HyperHDG/local_solver/diffusion_eigs_ldgh.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Configure eigenvector/-value solver.
  exact_eigenval = dimension * (np.pi ** 2)
  sigma          = exact_eigenval / 2.

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["shifted_inverse_eigenvalue_loop", problem, "vector[unsigned int]", \
            "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * 3 )
  helper = helper_ev_approx(HDG_wrapper, sigma)
  
  # Define LinearOperator in terms of C++ functions to use scipy's matrix-free linear solvers.
  system_size = helper.short_vector_size()
  Stiff       = LinearOperator( (system_size,system_size), matvec= helper.multiply_stiff )
  Mass        = LinearOperator( (system_size,system_size), matvec= helper.multiply_mass )
  ShiftedInv  = LinearOperator( (system_size,system_size), matvec= helper.shifted_inverse )
    
  # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
  vals, vecs = sp_lin_alg.eigsh(Stiff, k=1, M=Mass, sigma= sigma, which='LM', OPinv= ShiftedInv)

  # Print error.
  error = np.absolute(vals[0] - exact_eigenval)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_hypergraph_convergence_eigenvalue_shifted_inverse.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Postprocess solution vector.
  solution = helper.long_vector([x[0].real for x in vecs])
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "diff_eig_shifi-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(solution, vals[0].real)
  
  # Print ending time of diffusion test.
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  
  # Return smallest eigenvalue and corresponding eigenvector.
  return vals[0].real, solution, error
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(2,6):
        eigenvalue_approx_SI(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
