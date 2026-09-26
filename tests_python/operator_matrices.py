# Check trace_to_flux_mat (sparse_stiff_mat) and the operator cache (cache_operator) of all global
# loops that provide them, for the local solvers that implement trace_to_flux_mat themselves and for
# local solvers that use the fallback (trace_to_flux applied to unit vectors).
#
# For every case, the assembled matrix is compared to the matrix whose columns are the results of
# the matrix-free trace_to_flux applied to all unit vectors, and trace_to_flux with cached operator
# is compared to the matrix-free trace_to_flux.

import os, sys
import numpy as np
import scipy.sparse as sp

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG

def configure(global_loop, local_solver, parameters, dim, parabolic=False):
  const                 = HyperHDG.config()
  const.global_loop     = global_loop
  const.topology        = "Cubic<%d,%d>" % (dim, dim)
  const.geometry        = "UnitCube<%d,%d,double>" % (dim, dim)
  const.node_descriptor = "Cubic<%d,%d>" % (dim, dim)
  const.local_solver    = local_solver
  const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
  if parabolic:
    const.cython_replacements += ["double", "vector[double]"]
  const.include_files   = ["reproducibles_python/parameters/" + parameters]
  return HyperHDG.include(const)

# (name, global loop, local solver, parameter file, constructor arguments, parameter values,
#  test factorize_once (only for constant coefficients), parabolic)
CASES = [
  ("Bilaplacian", "Elliptic", "Bilaplacian<2,1,2,TestParametersSinEllipt,double>",
   "bilaplacian.hxx", {}, [0.], False, False),
  ("DiffusionIP", "Elliptic", "DiffusionIP<2,1,2,TestParametersSinEllipt,double>",
   "diffusion.hxx", {"lsol_constr": 20.}, [0.], True, False),
  ("DiffusionPostprocess", "Elliptic", "DiffusionPostprocess<2,2,4,TestParametersSinEllipt,double>",
   "diffusion.hxx", {}, [0.], True, False),
  ("DiffusionAdvectionReaction", "Elliptic",
   "DiffusionAdvectionReaction<2,1,2,TestParametersSinEllipt,double>",
   "diffusion_advection_reaction.hxx", {}, [0.], False, False),
  ("DiffusionParab", "Parabolic", "DiffusionParab<2,1,2,TestParametersSinParab,double>",
   "diffusion.hxx", {"lsol_constr": [1., 0.5, 0.1]}, [0., 0.3], True, True),
  ("BilaplacianParab", "Parabolic", "BilaplacianParab<2,1,2,TestParametersSinParab,double>",
   "bilaplacian.hxx", {"lsol_constr": [1., 1., 0.1]}, [0., 0.3], False, True),
  ("DiffusionEigs (fallback)", "ShiftedInverseEigenvalue",
   "DiffusionEigs<2,1,2,TestParametersEigs,double>", "diffusion.hxx", {}, [0., 5.], True, False),
  ("Diffusion", "MassApproxEigenvalue", "Diffusion<2,1,2,TestParametersSinEllipt,double>",
   "diffusion.hxx", {}, [0.], True, False),
]

# The elliptic global loop's Python interface has no time argument (the time is zero).
def apply(hdg, loop, vec, param):
  return np.asarray(hdg.trace_to_flux(vec) if loop == "Elliptic" else hdg.trace_to_flux(vec, param))

def assemble(hdg, loop, param, once):
  return hdg.sparse_stiff_mat(once) if loop == "Elliptic" else hdg.sparse_stiff_mat(param, once)

def dense_reference(hdg, loop, n, param):
  ref = np.zeros((n, n))
  for k in range(n):
    unit = np.zeros(n)
    unit[k] = 1.
    ref[:, k] = apply(hdg, loop, unit, param)
  return ref

success = True
for (name, loop, solver, params, constr, param_values, test_once, parabolic) in CASES:
  hdg = configure(loop, solver, params, 2, parabolic)([3, 2], **constr)
  n = hdg.size_of_system()
  x = np.random.default_rng(0).standard_normal(n)
  for param in param_values:
    ref = dense_reference(hdg, loop, n, param)
    scale = np.abs(ref).max()
    for once in ([False, True] if test_once else [False]):
      cols, rows, vals = assemble(hdg, loop, param, once)
      mat = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).toarray()
      dev = np.abs(mat - ref).max() / scale
      ok = scale > 0. and dev < 1e-12
      success &= ok
      print("%-28s %-25s param=%.2f factorize_once=%-5s deviation %.1e %s"
            % (name, loop, param, once, dev, "OK" if ok else "FAILED"))

  # Cached operator (also for changing parameter values, which require a new assembly).
  reference = [apply(hdg, loop, x, param) for param in param_values]
  hdg.cache_operator(True)
  cached = [apply(hdg, loop, x, param) for param in param_values * 2]
  hdg.cache_operator(False)
  dev = max(np.abs(c - r).max() / np.abs(r).max() for c, r in zip(cached, reference * 2))
  ok = dev < 1e-12
  success &= ok
  print("%-28s %-25s cache_operator                   deviation %.1e %s"
        % (name, loop, dev, "OK" if ok else "FAILED"))

assert success, "Some operator matrices differ from the matrix-free operators!"
