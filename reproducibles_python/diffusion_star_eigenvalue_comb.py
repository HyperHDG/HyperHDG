from __future__ import print_function

from diffusion_star_eigenvalue_approx import eigenvalue_approx_MA
from diffusion_star_eigenvalue_shifted_inverse import eigenvalue_approx_SI
from diffusion_star_eigenvalue_newt import eigenvalue_newt

import os, sys

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  for dimension in range(2, 3):
    print("\n Dimension is set to be ", dimension, "\n\n")
    for poly_degree in range(1, 3):
      print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
      for level in range(1, 4):
        print("\nLevel is ", level, "\n")
        for iteration in range(0, 3):
          for approx_type in range(0, 1):
            if approx_type == 0:
              value,vector,errorSI = eigenvalue_approx_SI(poly_degree,dimension,iteration,level,debug_mode)
              vector.append(value)
              try:
                _,_,errorNLSI = eigenvalue_newt(poly_degree,dimension,level,iteration, vector,debug_mode)
              except RuntimeError as error:
                print("ERROR: ", error)
                errorNLSI = -1.
              f = open("output/diffusion_star_eigenvalue_combSI.txt", "a")
              f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                      + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". ErrorSI = " + str(errorSI) + ".\n")
              f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                      + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLSI) + ".\n")
              f.close()
            else:
              value,vector,errorMA = eigenvalue_approx_MA(poly_degree,dimension,iteration,level,debug_mode)
              vector.append(value)
              try:
                _,_,errorNLMA = eigenvalue_newt(poly_degree,dimension,level,iteration,vector,debug_mode)
              except RuntimeError as error:
                print("ERROR: ", error)
                errorNLMA = -1.
              f = open("output/diffusion_star_eigenvalue_combMA.txt", "a")
              f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                      + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". ErrorMA = " + str(errorMA) + ".\n")
              f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
                      + ". Level = " + str(level) + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLMA) + ".\n")
              f.close()


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
