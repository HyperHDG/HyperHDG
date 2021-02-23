# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import other eigenvalue solvers.
from bilaplacian_hypergraph_convergence_eigenvalue_approx import eigenvalue_approx_MA
from bilaplacian_hypergraph_convergence_eigenvalue_shifted_inverse import eigenvalue_approx_SI
from bilaplacian_hypergraph_convergence_eigenvalue_newt import eigenvalue_newt

# Correct the python paths!
import os, sys
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for poly_degree in range(1,4):
    print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("\nDimension is ", dimension, "\n")
      for iteration in range(2,6):
        for approx_type in range(2):
          if approx_type == 0:
            value,vector,errorSI = eigenvalue_approx_SI(poly_degree,dimension,iteration,debug_mode)
            vector.append(value)
            try:
              _,_,errorNLSI = eigenvalue_newt(poly_degree,dimension,iteration,vector,debug_mode)
            except RuntimeError as error:
              print("ERROR: ", error)
              errorNLSI = -1.
            f = open("output/bilaplacian_hypergraph_convergence_eigenvalue_combSI.txt", "a")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorSI = " + str(errorSI) + ".\n")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLSI) + ".\n")
            f.close()
          else:
            value,vector,errorMA = eigenvalue_approx_MA(poly_degree,dimension,iteration,debug_mode)
            vector.append(value)
            try:
              _,_,errorNLMA = eigenvalue_newt(poly_degree,dimension,iteration,vector,debug_mode)
            except RuntimeError as error:
              print("ERROR: ", error)
              errorNLMA = -1.
            f = open("output/bilaplacian_hypergraph_convergence_eigenvalue_combMA.txt", "a")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorMA = " + str(errorMA) + ".\n")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLMA) + ".\n")
            f.close()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
