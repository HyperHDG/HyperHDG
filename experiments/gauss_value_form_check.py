#!/usr/bin/env python3
"""Order check of the eigenbasis Gauss step protocol used by timowave (HORK_GAUSS_PLAN.md).

Replicates, on a forced 2x2 oscillator u' = v + g1, v' = -u + g2 (exact solution
u = sin 2t, v = cos t), exactly the algorithm of the local solver / driver:

  per representative l:  (sigma_l I - J) zeta_l = sigma_l omega_l x^n + sum_j tinv_lj g(t_j)
  endpoint:              x+ = affine x^n + sum_l mult_l Re(w_l zeta_l)

with the tableau built like gauss_tableau.hxx (A = collocation, eigendecomposition,
omega = T^{-1} 1, w = (b^T T)/theta, affine = (-1)^s). Expected temporal order 2s for
s = 1, 2, 3 -- this isolates the tableau + recombination logic from all PDE machinery.
"""
import numpy as np

def build_tableau(s):
    c = (np.polynomial.legendre.leggauss(s)[0] + 1) / 2  # Gauss nodes on [0,1]
    P = np.vander(c, s, increasing=True)                 # P[j,k] = c_j^k
    rhs = np.array([[ci**(k+1) / (k+1) for k in range(s)] for ci in c])
    A = rhs @ np.linalg.inv(P)
    b = np.linalg.solve(P.T, np.array([1. / (k+1) for k in range(s)]))
    theta, T = np.linalg.eig(A.astype(complex))
    omega = np.linalg.solve(T, np.ones(s, dtype=complex))
    w = (b @ T) / theta
    affine = 1 - np.sum(w * omega)
    assert abs(affine.imag) < 1e-12 and abs(affine.real - (-1)**s) < 1e-10
    tinv = np.linalg.inv(T)
    # fold conjugate pairs: keep Im theta >= 0 representatives, real eigenvalues mult 1
    reps, mult = [], []
    for l in range(s):
        if theta[l].imag > 1e-12:
            reps.append(l); mult.append(2.)
        elif abs(theta[l].imag) <= 1e-12:
            reps.append(l); mult.append(1.)
    return c, theta, omega, w, tinv, affine.real, reps, np.array(mult)

def run(s, nt, T_end=1.):
    c, theta, omega, w, tinv, affine, reps, mult = build_tableau(s)
    J = np.array([[0., 1.], [-1., 0.]])
    u_ex = lambda t: np.sin(2*t)
    v_ex = lambda t: np.cos(t)
    g = lambda t: np.array([2*np.cos(2*t) - np.cos(t), -np.sin(t) + np.sin(2*t)])
    dt = T_end / nt
    x = np.array([u_ex(0.), v_ex(0.)])
    for n in range(nt):
        tn = n * dt
        loads = [g(tn + cj * dt) for cj in c]
        xp = affine * x.astype(complex)
        for l, m in zip(reps, mult):
            sigma = 1. / (theta[l] * dt)
            rhs = sigma * omega[l] * x + sum(tinv[l, j] * loads[j] for j in range(s))
            zeta = np.linalg.solve(sigma * np.eye(2, dtype=complex) - J, rhs)
            xp += m * (w[l] * zeta).real
        x = xp.real
    return np.hypot(x[0] - u_ex(T_end), x[1] - v_ex(T_end))

def run_prothero_robinson(s, nt, lam, T_end=1.):
    """Stiff order check: y' = lam*(y - phi(t)) + phi'(t), y(0) = phi(0), exact y = phi.

    Models the semidiscrete wave system's high spatial eigenmodes slaved to smooth
    time-dependent boundary data (wave4): for |lam*dt| >~ 1 Gauss collocation loses its
    endpoint superconvergence (classical order reduction, expect ~s+1 instead of 2s).
    """
    c, theta, omega, w, tinv, affine, reps, mult = build_tableau(s)
    phi = lambda t: np.cos(2*t + 0.3)
    dphi = lambda t: -2*np.sin(2*t + 0.3)
    g = lambda t: -lam * phi(t) + dphi(t)
    dt = T_end / nt
    y = phi(0.) + 0j
    for n in range(nt):
        tn = n * dt
        yp = affine * y
        for l, m in zip(reps, mult):
            sigma = 1. / (theta[l] * dt)
            rhs = sigma * omega[l] * y + sum(tinv[l, j] * g(tn + c[j] * dt) for j in range(s))
            zeta = rhs / (sigma - lam)
            yp += m * (w[l] * zeta).real
        y = yp.real
    return abs(y - phi(T_end))

for s in (1, 2, 3):
    errs = [run(s, nt) for nt in (4, 8, 16, 32, 64)]
    orders = [np.log2(errs[k-1] / errs[k]) for k in range(1, len(errs))]
    print(f"s={s}: errs={['%.3e' % e for e in errs]} orders={['%.2f' % o for o in orders]}")

for lam in (1e3j, -1e3, 1e6j):
    for s in (2, 3):
        errs = [run_prothero_robinson(s, nt, lam) for nt in (4, 8, 16, 32, 64)]
        orders = [np.log2(errs[k-1] / errs[k]) for k in range(1, len(errs))]
        print(f"PR lam={lam}: s={s}: errs={['%.3e' % e for e in errs]} "
              f"orders={['%.2f' % o for o in orders]}")
