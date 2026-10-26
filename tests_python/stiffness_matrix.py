# Check that sparse_stiff_mat yields the same matrix with and without factorize_once on a uniform
# mesh with constant coefficients, and that the resulting system can be solved.

import os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG

dimension, poly_degree, n_elements = 2, 2, 8

const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
const.geometry        = "UnitCube<" + str(dimension) + "," + str(dimension) + ",double>"
const.node_descriptor = "Cubic<" + str(dimension) + "," + str(dimension) + ">"
const.local_solver    = "Diffusion<" + str(dimension) + "," + str(poly_degree) + "," \
  + str(2*poly_degree) + ",TestParametersSinEllipt,double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
const.debug_mode      = False

HDG_wrapper = HyperHDG.include(const)( [n_elements] * dimension )
n = HDG_wrapper.size_of_system()

matrices = []
for args in [(), (False,), (True,)]:
  col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(*args)
  # The triplets are returned as NumPy arrays of matching lengths.
  assert all(isinstance(array, np.ndarray) for array in (col_ind, row_ind, vals))
  assert col_ind.dtype == row_ind.dtype == np.uintc and vals.dtype == np.float64
  assert len(col_ind) == len(row_ind) == len(vals) > 0
  matrices.append(sp.csr_matrix((vals, (row_ind, col_ind)), shape=(n, n)))

for A in matrices[1:]:
  assert abs(A - matrices[0]).max() <= 1e-14 * abs(matrices[0]).max(), "Matrices differ!"

rhs = np.multiply(HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1.)
solution, info = sp_lin_alg.cg(matrices[2], rhs, tol=1e-13)
assert info == 0, "CG did not converge!"
error = HDG_wrapper.errors(solution)[0]
print("Error with factorize_once = True:", error)
assert error < 1e-3, "Error is too large!"

# On a non-uniform hypergraph, factorize_once is wrong. This must raise a ValueError (the Python
# interface is compiled in release mode, where only a sample of hyperedges is checked). In the
# refined triangle, the interior hyperedges of the legs and of the hypotenuse have the same face
# types, but different lengths.
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "File<1,2>"
const.geometry        = "File<1,2>"
const.node_descriptor = "File<1,2>"
const.local_solver    = "Diffusion<1,2,4,TestParametersSinEllipt,double>"
const.cython_replacements = ["string", "string"]
const.include_files   = ["reproducibles_python/parameters/diffusion.hxx"]
const.debug_mode      = False

HDG_wrapper = HyperHDG.include(const)(os.path.dirname(os.path.abspath(__file__)) \
  + "/../domains/simplex_1_2.geo")
HDG_wrapper.refine(4)
HDG_wrapper.sparse_stiff_mat(False)
try:
  HDG_wrapper.sparse_stiff_mat(True)
  raise AssertionError("Misuse of factorize_once has not been detected!")
except ValueError as error:
  print("Misuse of factorize_once detected:", str(error)[:60], "...")
