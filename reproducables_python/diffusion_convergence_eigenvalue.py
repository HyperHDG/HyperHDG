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
class helper_class():
  def __init__(self, hdg_wrapper):
    self.hdg_wrapper = hdg_wrapper
    self.val = [0] * (hdg_wrapper.size_of_system() + 1)
  def eval_residual(self, vector):
    vec = self.hdg_wrapper.matrix_vector_multiply(vector, vector[len(vector)-1])
    temp = vector[len(vector)-1]
    vector[len(vector)-1] = 0.
    vec[len(vec)-1] = np.linalg.norm(vector) ** 2 - 1.
    vector[len(vector)-1] = temp
    return vec
  def set_val(self, val):
    self.val = val;
  def eval_jacobian(self, vector):
    vec = self.hdg_wrapper.matrix_vector_der_multiply \
            ( vector, vector[len(vector)-1], self.val, self.val[len(self.val)-1] )
    temp = self.val[len(self.val)-1]
    self.val[len(self.val)-1] = 0.
    vec[len(vec)-1] = 2. * np.dot(vector,self.val)
    self.val[len(self.val)-1] = temp
    return vec


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def diffusion_test(poly_degree, dimension, iteration):
  
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "DiffusionEigs<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersHomo,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/Diffusion.hxx", \
                "reproducables_python/parameters/diffusion.hxx" ]

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )
  helper = helper_class(HDG_wrapper)

  # Define LinearOperator in terms of C++ functions to use scipy in a matrix-free fashion.
  system_size = HDG_wrapper.size_of_system() + 1
  A = LinearOperator( (system_size,system_size), matvec= helper.eval_jacobian )
  
  # Initialize solution vector [lambda, eig].
  vectorSolution = [0] * system_size
  
  # Initial vector is solution!
  for i in range (system_size-2):
    vectorSolution[i] = np.sin(np.pi * i / (system_size - 2))
  vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
  
  # Config Newton solver
  n_newton_steps = 10 ** 1
  norm_res       = 1e-9
  vectorSolution[system_size-1] = np.pi * np.pi
  
  # residual = helper.eval_residual(vectorSolution)
  
  # print(np.linalg.norm( residual ))
  
  # For loop over the respective time-steps.
  for newton_step in range(n_newton_steps):
    
    #helper.set_val(vectorSolution)
    
    # print(residual)
    # print(A * residual) 
    
    # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
    # [vectorUpdate, num_iter] = sp_lin_alg.cg(A, residual, tol=1e-13)
    # if num_iter != 0:
    #   print("CG failed with a total number of ", num_iter, " iterations in time step ", newton_step, \
    #        ". Trying GMRES!")
    #  [vectorUpdate, num_iter] = sp_lin_alg.gmres(A, residual, tol=1e-13)
    #  if num_iter != 0:
    
    #    print("GMRES also failed with a total number of ", num_iter, "iterations.")
    #    sys.exit("Program failed!")
    
    # print(vectorUpdate)
    # vectorSolution = np.subtract(vectorSolution, vectorUpdate)
    
    residual = helper.eval_residual(vectorSolution)
    #print(residual)
    vectorSolution = np.add(vectorSolution, np.multiply(residual,1.) )
    
    # print(residual)
    print(np.linalg.norm( residual ))
    
    # residual = helper.eval_residual(vectorSolution)
    # if np.linalg.norm( residual ) < norm_res:
    #  break

  # Print error.
  error = np.absolute(vectorSolution[system_size-1] - np.pi * np.pi)
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/diffusion_convergence_eigenvalue.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  #HDG_wrapper.plot_option( "fileName" , "diff_e-" + str(dimension) + "-" + str(iteration) );
  #HDG_wrapper.plot_option( "printFileNumber" , "false" );
  #HDG_wrapper.plot_option( "scale" , "0.95" );
  #HDG_wrapper.plot_solution(HDG_wrapper.plot_solution(vectorSolution,vectorSolution[system_size-1]))
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,2):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,2):
      print("Dimension is ", dimension, "\n")
      for iteration in range(5,6):
        diffusion_test(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
