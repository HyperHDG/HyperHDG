# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import other eigenvalue solvers.
from diffusion_convergence_eigenvalue_approx import eigenvalue_approx
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
        value, vector = eigenvalue_approx(poly_degree, dimension, iteration)
        vector.append(value)
        eigenvalue_newt(poly_degree, dimension, iteration, vector)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
