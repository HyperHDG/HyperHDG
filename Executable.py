# Python running example

from __future__ import print_function

import numpy as np
import scipy.sparse.linalg as sp_lin_alg

from ClassWrapper import PyDiffusionProblem
from scipy.sparse.linalg import LinearOperator

# Initialising the wrapped c++ function
R1 = PyDiffusionProblem([2,2], 1);

vectorDirichlet = R1.return_zero_vector();
vectorDirichlet[0] = 1.;
vectorDirichlet[8] = 0.;

print("Test 1: ", vectorDirichlet)

vectorMultiply = R1.matrix_vector_multiply(vectorDirichlet);
vectorMultiply = [-i for i in vectorMultiply];
print("Test 2: ", vectorMultiply)

system_size = R1.size_of_system();
A = LinearOperator((system_size,system_size), matvec= R1.matrix_vector_multiply);

vector_x = sp_lin_alg.cg(A, vectorMultiply, maxiter=100);
print("Test 4: ", vector_x)
