from __future__ import print_function

from timoshenko_convergence_eigenvalue_shifted_inverse import eigenvalue_approx_SI
from timoshenko_convergence_eigenvalue_newt import eigenvalue_newt

import os, sys
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  dimension = 1
  for poly_degree in range(4,5):
    print("\n  Polynomial degree is set to be ", poly_degree, "\n\n")
    for iteration in range(1,7):
      value,vector,errorSI = eigenvalue_approx_SI(poly_degree,dimension,iteration,debug_mode)
      vector.append(value)
      try:
        _,_,errorNLSI = eigenvalue_newt(poly_degree,dimension,iteration,vector,debug_mode)
      except RuntimeError as error:
        print("ERROR: ", error)
        errorNLSI = -1.
      f = open("output/diffusion_hypergraph_convergence_eigenvalue_combSI.txt", "a")
      f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
              + ". Iteration = " + str(iteration) + ". ErrorSI = " + str(errorSI) + ".\n")
      f.write("Polynomial degree = " + str(poly_degree) + ". Dimension = " + str(dimension) \
              + ". Iteration = " + str(iteration) + ". ErrorNL = " + str(errorNLSI) + ".\n")
      f.close()

# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
