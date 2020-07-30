# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import numpy and spares linear algebra from scipy to use the Python maths libraries.
import numpy as np
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator
import random

# Correct the python paths!
import os, sys
sys.path.append(os.path.dirname(__file__) + "/..")


# --------------------------------------------------------------------------------------------------
# Class implementing the matvec of "mass_matrix + delta_time * stiffness_matrix".
# --------------------------------------------------------------------------------------------------
class helper_ev_newt():
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
def eigenvalue_newt(poly_degree, dimension, iteration, initial="default"):
  
  # Predefine problem to be solved.
  problem = "AbstractProblem < Topology::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "Geometry::UnitCube<" + str(dimension) + "," + str(dimension) + ",double>, " \
          + "NodeDescriptor::Cubic<" + str(dimension) + "," + str(dimension) + ">, " \
          + "bilaplacian_eigs<" + str(dimension) + "," + str(poly_degree) + "," + str(2*poly_degree) \
          + ",TestParametersEigs,double> >"
  filenames = [ "HyperHDG/Geometry/Cubic.hxx" , "HyperHDG/NodeDescriptor/Cubic.hxx", \
                "HyperHDG/LocalSolver/bilaplacian.hxx", \
                "reproducables_python/parameters/bilaplacian.hxx" ]

  # Config Newton solver
  n_newton_steps = 10 ** 1
  scaling_fac    = 1e-2
  norm_exact     = 1e-9
  alpha          = 0.1
  beta           = 0.5

  # Import C++ wrapper class to use HDG method on graphs.
  from cython_import import cython_import
  PyDP = cython_import \
         ( ["AbstractProblem", problem, "vector[unsigned int]", "vector[unsigned int]"], filenames )

  # Initialising the wrapped C++ class HDG_wrapper.
  HDG_wrapper = PyDP( [2 ** iteration] * dimension )
  helper = helper_ev_newt(HDG_wrapper)

  # Define LinearOperator in terms of C++ functions to use scipy linear solvers in a matrix-free
  # fashion.
  system_size = HDG_wrapper.size_of_system() + 1
  A = LinearOperator( (system_size,system_size), matvec= helper.eval_jacobian )
  
  # Initialize solution vector [lambda, eig].
  if initial == "default":
    vectorSolution = [0] * system_size
    vectorSolution = HDG_wrapper.initial_flux_vector(vectorSolution)
    vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    vectorSolution[system_size-1] = dimension * (np.pi ** 4) + 1e-3 * random.randint(-100,100)
  else:
    vectorSolution = initial
    temp = initial[len(vectorSolution)-1]
    vectorSolution[len(vectorSolution)-1] = 0.
    vectorSolution = np.multiply(vectorSolution, 1./np.linalg.norm(vectorSolution))
    vectorSolution[len(vectorSolution)-1] = temp
  
  # Initial residual.
  residual = helper.eval_residual( vectorSolution )
  norm_res = np.linalg.norm( residual )
  
  # For loop over the respective time-steps.
  for newton_step in range(n_newton_steps):
    
    helper.set_val(vectorSolution)
    norm_old = norm_res
    gamma = 1
    
    # Solve "A * x = b" in matrix-free fashion using scipy's CG algorithm.
    # [vectorUpdate, num_iter] = sp_lin_alg.cg(A, residual, tol=1e-13)
    # if num_iter != 0:
    #   print("CG failed with a total number of ", num_iter, " iterations in time step ", \
    #         newton_step, ". Trying GMRES!")
    [vectorUpdate, num_iter] = sp_lin_alg.bicgstab(A,residual,tol=min(1e-9,scaling_fac * norm_res))
    if num_iter != 0:
      print("BiCGStab failed with a total number of ", num_iter, "iterations.")
      [vectorUpdate, num_iter] = sp_lin_alg.gmres(A,residual,tol=min(1e-9,scaling_fac * norm_res))
      if num_iter != 0:
        print("GMRES also failed with a total number of ", num_iter, "iterations.")
        sys.exit("Program failed!")
    
    vectorHelper = np.subtract(vectorSolution, vectorUpdate)
    residual = helper.eval_residual(vectorHelper)
    norm_res = np.linalg.norm( residual )
    
    while norm_res > (1 - alpha * gamma) * norm_old:
      if gamma < 1e-4:
        sys.exit("Newton step is too small!")
      gamma = beta * gamma
      vectorHelper = np.subtract(vectorSolution, np.multiply(vectorUpdate, gamma))
      residual = helper.eval_residual(vectorHelper)
      norm_res = np.linalg.norm( residual )
    
    vectorSolution = vectorHelper
    
    if norm_res < norm_exact:
      print("Newton solver converged after ", newton_step+1, " steps with a residual of norm ",\
            norm_res)
      break  
  
  if norm_res >= norm_exact:
    print("Newton solver did not converge with final residual = ", norm_res)
    sys.exit("Program failed!")

  # Print error.
  error = np.absolute(vectorSolution[system_size-1] - dimension * (np.pi ** 4))
  print("Iteration: ", iteration, " Error: ", error)
  f = open("output/bilaplacian_convergence_eigenvalue.txt", "a")
  f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
          + ". Iteration = " + str(iteration) + ". Error = " + str(error) + ".\n")
  f.close()
  
  # Plot obtained solution.
  HDG_wrapper.plot_option( "fileName" , "bil_e-" + str(dimension) + "-" + str(iteration) )
  HDG_wrapper.plot_option( "printFileNumber" , "false" )
  HDG_wrapper.plot_option( "scale" , "0.95" )
  HDG_wrapper.plot_solution(vectorSolution,vectorSolution[system_size-1])
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("Dimension is ", dimension, "\n")
      for iteration in range(1,6):
        bilaplacian_test(poly_degree, dimension, iteration)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()