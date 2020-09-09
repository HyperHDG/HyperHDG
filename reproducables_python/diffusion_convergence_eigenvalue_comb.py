# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import other eigenvalue solvers.
from diffusion_convergence_eigenvalue_approx import eigenvalue_approx_MA
from diffusion_convergence_eigenvalue_shifted_inverse import eigenvalue_approx_SI
from diffusion_convergence_eigenvalue_newt import eigenvalue_newt

# Correct the python paths!
import os, sys
sys.path.append(os.path.dirname(__file__) + "/..")
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main():
  for poly_degree in range(1,4):
    print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
    for dimension in range(1,3):
      print("\nDimension is ", dimension, "\n")
      for iteration in range(2,6):
        for approx_type in range(2):
          if approx_type == 0:
            value, vector, errorSI = eigenvalue_approx_SI(poly_degree, dimension, iteration)
            vector.append(value)
            _, _, errorNLSI = eigenvalue_newt(poly_degree, dimension, iteration, vector)
            f = open("output/diffusion_convergence_eigenvalue_combSI.txt", "a")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorSI = " + str(errorSI) + ".\n")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLSI) + ".\n")
            f.close()
          else:
            value, vector, errorMA = eigenvalue_approx_MA(poly_degree, dimension, iteration)
            vector.append(value)
            _, _, errorNLMA = eigenvalue_newt(poly_degree, dimension, iteration, vector)
            f = open("output/diffusion_convergence_eigenvalue_combMA.txt", "a")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorMA = " + str(errorMA) + ".\n")
            f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                    + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLMA) + ".\n")
            f.close()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
