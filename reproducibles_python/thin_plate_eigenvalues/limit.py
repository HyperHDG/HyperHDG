"""Exact eigenpairs of the limit problem of three unit-square plates meeting along a straight line.

Plate E_i is parametrized by s (along the junction N, 0 < s < 1) and r (distance from N,
0 < r < 1). The plate weights are d = (3/2, 3/2, 3), kappa_i = rho_i = 1, the junction carries the
mass density mu = 24 / (17 pi), the plates are clamped (Dirichlet) at r = 1, and the edges s = 0
and s = 1 of all plates are natural boundaries. Separation of variables gives the eigenfunctions
u^i(s, r) = f^i(r) cos(k pi s), k = 0, 1, 2, ..., where either
  - (symmetric) f^i(r) = sin(omega (1 - r)) on all plates with
    sum(d) omega cos(omega) = (omega^2 + k^2 pi^2) mu sin(omega), lambda = omega^2 + k^2 pi^2, or
  - (antisymmetric) f^i(r) = a_i sin(n pi r) with d_1 a_1 + d_2 a_2 + d_3 a_3 = 0, u^N = 0, and
    lambda = (n^2 + k^2) pi^2 (double eigenvalue).
"""

import math

import numpy as np
from scipy.optimize import brentq

D = np.array([1.5, 1.5, 3.0])
MU = 24. / (17. * math.pi)


def symmetric_omega(k, branch=0, d=D, mu=MU):
  """Root number branch (0, 1, ...) of sum(d) w cos w = (w^2 + k^2 pi^2) mu sin w, w > 0."""
  func = lambda w: d.sum() * w * math.cos(w) - (w * w + (k * math.pi) ** 2) * mu * math.sin(w)
  grid = np.linspace(1e-9, (branch + 2) * math.pi, 20000 * (branch + 2) + 1)
  values = [func(w) for w in grid]
  roots = [brentq(func, a, b, xtol=1e-15, rtol=1e-15) for a, b, fa, fb in
           zip(grid[:-1], grid[1:], values[:-1], values[1:]) if fa * fb < 0.]
  roots = [w for w in roots if abs(math.sin(w)) > 1e-12]
  return roots[branch]


class Mode:
  """Eigenfunction u^i(s, r) = a_i g(r) cos(k pi s) of the limit problem (with junction value)."""

  def __init__(self, k, omega, coeffs, junction_value):
    self.k, self.omega = k, omega
    self.coeffs = np.asarray(coeffs, dtype=float)
    self.junction_value = junction_value  # u^N(s) = junction_value * cos(k pi s)
    self.eigenvalue = omega ** 2 + (k * math.pi) ** 2

  def radial(self, plate, r):
    return self.coeffs[plate] * np.sin(self.omega * (1. - np.asarray(r)))

  def axial(self, s):
    return np.cos(self.k * math.pi * np.asarray(s))

  def __call__(self, plate, s, r):
    return self.radial(plate, r) * self.axial(s)


def symmetric_mode(k, branch=0):
  omega = symmetric_omega(k, branch)
  return Mode(k, omega, [1., 1., 1.], math.sin(omega))


def antisymmetric_mode(k, n, coeffs):
  # sin(n pi r) = (-1)^(n+1) sin(n pi (1 - r)); the sign is irrelevant for eigenspaces.
  return Mode(k, n * math.pi, coeffs, 0.)


def clusters():
  """The three clusters studied: lambda_1 (simple), lambda_2 = lambda_3 = pi^2, lambda_4."""
  return {
    "lambda1": {"first_index": 0, "modes": [symmetric_mode(0)]},
    "lambda23": {"first_index": 1, "modes": [antisymmetric_mode(0, 1, [1., -1., 0.]),
                                             antisymmetric_mode(0, 1, [1., 1., -1.])]},
    "lambda4": {"first_index": 3, "modes": [symmetric_mode(1)]},
  }


def eigenvalue(cluster):
  values = [mode.eigenvalue for mode in cluster["modes"]]
  assert max(values) - min(values) < 1e-12
  return values[0]


def radial_inner(mode1, mode2, n_quad=60):
  """sum_i d_i int_0^1 f_1^i f_2^i dr + mu f_1^N f_2^N (the s-factor is not included)."""
  points, weights = np.polynomial.legendre.leggauss(n_quad)
  r, w = 0.5 * (points + 1.), 0.5 * weights
  value = sum(D[i] * np.dot(w, mode1.radial(i, r) * mode2.radial(i, r)) for i in range(3))
  return value + MU * mode1.junction_value * mode2.junction_value


def axial_inner(mode1, mode2):
  """int_0^1 cos(k1 pi s) cos(k2 pi s) ds."""
  if mode1.k != mode2.k:
    return 0.
  return 1. if mode1.k == 0 else 0.5


def h0_inner(mode1, mode2):
  """Inner product of H_0 (plate and junction mass) of two limit eigenfunctions."""
  return radial_inner(mode1, mode2) * axial_inner(mode1, mode2)


def gram(modes, inner=h0_inner):
  return np.array([[inner(a, b) for b in modes] for a in modes])


def spectrum(n_eigs=8):
  """The n_eigs smallest eigenvalues (with multiplicity) and a description of their modes."""
  entries = []
  for k in range(4):
    for branch in range(3):
      entries.append((symmetric_omega(k, branch) ** 2 + (k * math.pi) ** 2, f"k={k} symmetric"))
    for n in range(1, 4):
      entries += 2 * [((n * n + k * k) * math.pi ** 2, f"k={k} antisymmetric n={n}")]
  # imaginary frequencies (lambda < k^2 pi^2) only exist for k >= 2, i.e., above 4 pi^2 - 6/mu
  return sorted(entries)[:n_eigs]


if __name__ == "__main__":
  for index, (value, kind) in enumerate(spectrum(8)):
    print(f"lambda_{index + 1} = {value:.15f}  ({kind})")
