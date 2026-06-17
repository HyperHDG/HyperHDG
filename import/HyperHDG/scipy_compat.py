## \package HyperHDG.scipy_compat
#
#  \brief   Keep the legacy ``tol`` keyword of SciPy's iterative solvers working.
#
#  SciPy 1.12 renamed the relative-tolerance keyword of its sparse iterative solvers (``cg``,
#  ``gmres``, ``bicgstab``, ...) from ``tol`` to ``rtol`` and removed ``tol`` entirely in SciPy
#  1.14. The HyperHDG example and reproducer scripts were written against the old name. Importing
#  the HyperHDG package installs the thin wrappers below so that those scripts keep working on
#  recent SciPy while remaining untouched on older versions.
#
#  The patch is intentionally conservative: a solver is only wrapped when the installed SciPy no
#  longer accepts ``tol`` but does accept ``rtol``. On SciPy < 1.12 (e.g. the versions shipped by
#  some Linux distributions) every solver still has a native ``tol`` parameter, so nothing is
#  wrapped and behaviour is unchanged.

import functools
import inspect

## Iterative solvers that received the tol -> rtol rename.
_AFFECTED_SOLVERS = ("bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres", "minres", "qmr",
                     "gcrotmk", "tfqmr")

## \brief   Wrap a solver so that a passed ``tol`` keyword is forwarded as ``rtol``.
def _tol_to_rtol(func):
  try:
    params = inspect.signature(func).parameters
  except (TypeError, ValueError):
    return func
  # Only wrap when the modern signature dropped ``tol`` in favour of ``rtol``.
  if "rtol" not in params or "tol" in params:
    return func

  @functools.wraps(func)
  def wrapper(*args, **kwargs):
    if "tol" in kwargs and "rtol" not in kwargs:
      kwargs["rtol"] = kwargs.pop("tol")
    return func(*args, **kwargs)
  return wrapper

## \brief   Patch ``scipy.sparse.linalg`` in place so the legacy ``tol`` keyword keeps working.
def patch():
  try:
    import scipy.sparse.linalg as ssl
  except (ImportError, ModuleNotFoundError):
    return
  for name in _AFFECTED_SOLVERS:
    func = getattr(ssl, name, None)
    if func is not None:
      setattr(ssl, name, _tol_to_rtol(func))
