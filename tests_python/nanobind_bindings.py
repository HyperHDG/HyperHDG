# Check the explicit nanobind bindings (include/HyperHDG/bind_python.hxx, python/hyperhdg.py).
#
# - The module is built for the running interpreter, and the built extension module (.so or .pyd)
#   is found and imported.
# - trace_to_flux_mat returns NumPy arrays and agrees with the matrix-free trace_to_flux, also with
#   factorize_once.
# - The diffusion problem is solved by CG with the assembled matrix.
# - Two variants of a module get separate build directories, so switching between them does not
#   rebuild them.

import inspect, os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../python")
import hyperhdg

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

using HDG = GlobalLoop::Elliptic<
  Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>, NodeDescriptor::Cubic<2, 2>,
  LocalSolver::Diffusion<2, {degree}, {quadrature}, TestParametersSinEllipt, double> >;

NB_MODULE(hdg_test_diffusion, m) {{
  HyperHDG::bind_python<HDG, std::vector<unsigned int> >(m, "Diffusion");
}}
"""

# SciPy renamed the tolerance of cg from tol to rtol in version 1.12 (Ubuntu 24.04 has 1.11).
TOLERANCE = "rtol" if "rtol" in inspect.signature(sp_lin_alg.cg).parameters else "tol"

def code(degree):
  return MODULE.format(degree=degree, quadrature=2 * degree)

def solve(hdg):
  n = hdg.size_of_system()
  rows, cols, vals = hdg.trace_to_flux_mat()
  A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
  rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
  solution, info = sp_lin_alg.cg(A, rhs, **{TOLERANCE: 1e-13})
  assert info == 0, "CG did not converge!"
  return hdg.errors(list(solution))[0]

success = True
module = hyperhdg.load(code(1))
hdg = module.Diffusion([4, 4], 1.)
n = hdg.size_of_system()

# The triplets are NumPy arrays and represent the matrix-free operator.
reference = np.zeros((n, n))
for k in range(n):
  unit = np.zeros(n)
  unit[k] = 1.
  reference[:, k] = hdg.trace_to_flux(list(unit))
scale = np.abs(reference).max()
for once in [False, True]:
  rows, cols, vals = hdg.trace_to_flux_mat(0., once)
  is_numpy = all(isinstance(array, np.ndarray) for array in (rows, cols, vals))
  dense = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).toarray()
  deviation = np.abs(dense - reference).max() / scale
  ok = is_numpy and scale > 0. and deviation < 1e-12
  success &= ok
  print("trace_to_flux_mat factorize_once=%-5s NumPy %s deviation %.1e %s"
        % (once, is_numpy, deviation, "OK" if ok else "FAILED"))

# The solution converges (the error is reduced by refinement).
errors = [solve(module.Diffusion([2 ** level] * 2, 1.)) for level in [2, 3]]
ok = errors[1] < 0.5 * errors[0] < 1e-1
success &= ok
print("errors", errors, "OK" if ok else "FAILED")

# Two variants are built in separate directories, and both stay valid.
paths = [hyperhdg.compile(code(degree), hash_name=True) for degree in [1, 2, 1]]
ok = paths[0] == paths[2] and paths[0].parent != paths[1].parent and paths[0].is_file() \
     and paths[1].is_file()
success &= ok
print("separate build directories for variants:", "OK" if ok else "FAILED")

assert success, "Some checks of the nanobind bindings failed!"
