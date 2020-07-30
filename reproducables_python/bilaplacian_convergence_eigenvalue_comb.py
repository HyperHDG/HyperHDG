# Python running example for HDG solution of an elasticity problem on a superaggregate!

# Import printing functions.
from __future__ import print_function

# Import other eigenvalue solvers.
from bilaplacian_convergence_eigenvalue_approx import eigenvalue_approx
from bilaplacian_convergence_eigenvalue_newt import eigenvalue_newt

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
        initial = [0.] * (len(vector) + 1)
        for i in range(len(vector)):
          initial[i] = (vector[i]).real
        initial[len(vector)] = value.real
        eigenvalue_newt(poly_degree, dimension, iteration, initial)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
    main()
