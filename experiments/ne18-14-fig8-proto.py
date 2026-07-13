#!/usr/bin/env python
"""Verbatim scipy prototype of gortz.pdf's subspace-decomposition PCG (sec. 4) for the
regular-grid heat problem (sec. 6.2, Table 2 row "Grid"), to isolate implementation
differences from net2as cb_q1: interior dofs only (V with u=0 on the boundary), coarse
space Q_H = BC-conforming interior Q1 hats, one subspace V(U(y_j)) per coarse mesh node
(all nodes, closed element patches U), exact patch/coarse solves, plain additive PCG.

Prints per-iteration K-norm errors, (tau_avg, tau_max) as in eq. (6.2), and the extreme
eigenvalues of BK estimated from the CG Lanczos tridiagonal (as PETSc does)."""

import argparse

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=513, help="grid nodes per axis (2^9+1)")
    ap.add_argument("--Hinv", type=int, default=8, help="coarse cells per axis, H = 1/Hinv")
    ap.add_argument("--iters", type=int, default=60)
    ap.add_argument("--open-patches", action="store_true",
                    help="open hat supports instead of closed element patches U(y_j)")
    ap.add_argument("--keep-bnd-cb", action="store_true",
                    help="keep the boundary hats in the coarse space (net2as default)")
    args = ap.parse_args()

    n, Hn = args.n, args.Hinv
    h, H = 1.0 / (n - 1), 1.0 / args.Hinv
    m = n - 2                       # interior nodes per axis
    N = m * m

    # K = weighted graph Laplacian (conductance 1/h), Dirichlet boundary eliminated
    T = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(m, m))
    I = sp.identity(m, format="csr")
    K = ((sp.kron(I, T) + sp.kron(T, I)) / h).tocsr()

    x1 = np.arange(1, n - 1) * h
    X, Y = np.meshgrid(x1, x1, indexing="ij")
    xs, ys = X.ravel(), Y.ravel()

    def hat(a, b):
        wx = np.clip(1.0 - np.abs(xs - a * H) / H, 0.0, None)
        wy = np.clip(1.0 - np.abs(ys - b * H) / H, 0.0, None)
        return wx * wy

    # coarse basis
    cols = []
    for a in range(Hn + 1):
        for b in range(Hn + 1):
            interior = 0 < a < Hn and 0 < b < Hn
            if interior or args.keep_bnd_cb:
                cols.append(hat(a, b))
    R0 = sp.csr_matrix(np.array(cols))
    A0inv = np.linalg.inv((R0 @ K @ R0.T).toarray())

    # subdomain patches: one per coarse mesh node (all of them)
    patches, lus = [], []
    for a in range(Hn + 1):
        for b in range(Hn + 1):
            if args.open_patches:
                msk = (np.abs(xs - a * H) < H) & (np.abs(ys - b * H) < H)
            else:  # closed element patch U(y_j): elements whose closure touches y_j
                tol = 1e-12
                msk = (np.abs(xs - a * H) <= H + tol) & (np.abs(ys - b * H) <= H + tol)
            ids = np.nonzero(msk)[0]
            if len(ids):
                patches.append(ids)
                lus.append(spla.splu(K[ids][:, ids].tocsc()))
    print(f"n={n} H=1/{Hn}: {N} interior dofs, {R0.shape[0]} coarse dofs, "
          f"{len(patches)} patches (sizes {min(map(len, patches))}..{max(map(len, patches))})")

    def B(r):
        y = R0.T @ (A0inv @ (R0 @ r))
        for ids, lu in zip(patches, lus):
            y[ids] += lu.solve(r[ids])
        return y

    # heat problem: f = M1 (node mass = half the length of adjacent edges)
    b = np.full(N, 2.0 * h)
    uref = spla.spsolve(K.tocsc(), b)

    # PCG in the K-inner product, recording K-norm errors and Lanczos coefficients
    x = np.zeros(N)
    r = b.copy()
    z = B(r)
    p = z.copy()
    rz = r @ z
    e = uref - x
    errs = [np.sqrt(e @ (K @ e))]
    alphas, betas = [], []
    for _ in range(args.iters):
        Kp = K @ p
        alpha = rz / (p @ Kp)
        x += alpha * p
        r -= alpha * Kp
        e = uref - x
        errs.append(np.sqrt(max(e @ (K @ e), 0.0)))
        alphas.append(alpha)
        if errs[-1] < 1e-10 * errs[0]:
            break
        z = B(r)
        rz_new = r @ z
        betas.append(rz_new / rz)
        p = z + (rz_new / rz) * p
        rz = rz_new

    errs = np.array(errs)
    rates = errs[1:] / errs[:-1]
    tau_avg = np.exp(np.log(errs[-1] / errs[1]) / (len(errs) - 2))
    print("errs:", " ".join(f"{v:.3e}" for v in errs))
    print(f"iters: {len(errs) - 1}  tau_avg: {tau_avg:.3f}  tau_max: {rates[1:].max():.3f}")

    # extreme eigenvalues of BK from the CG Lanczos tridiagonal
    k = len(alphas)
    diag = np.array([1.0 / alphas[0]] +
                    [1.0 / alphas[i] + betas[i - 1] / alphas[i - 1] for i in range(1, k)])
    off = np.array([np.sqrt(betas[i]) / alphas[i] for i in range(k - 1)])
    ev = np.linalg.eigvalsh(np.diag(diag) + np.diag(off, 1) + np.diag(off, -1))
    print(f"emin: {ev[0]:.4f}  emax: {ev[-1]:.4f}  cond: {ev[-1] / ev[0]:.3f}")


if __name__ == "__main__":
    main()
