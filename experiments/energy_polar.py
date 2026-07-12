#!/usr/bin/env python3
"""Angular energy distribution E(theta, t) from a wave.vtkhdf with per-cell energies.

Companion to energy_rt.py (same input, same energy definition): bins the per-cell
energy over the polar angle of the cell center around the tap center, and plots
  (1) the E(theta, t) heatmap -- lattice anisotropy shows as vertical banding
      (a square grid runs faster along its axes than its diagonals), and
  (2) polar snapshots E(theta) at selected times, restricted to an annulus around the
      instantaneous front radius so the anisotropy of the FRONT is not washed out by
      the isotropic source region.
"""

import argparse
import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from energy_rt import read_cell_energies


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file", help="path to wave.vtkhdf (needs CellData/energies)")
    ap.add_argument("-o", "--save", default="energy_polar.png", help="output image path")
    ap.add_argument("--center", default=None, metavar="X,Y",
                    help="tap center (default: xy bbox center of the mesh)")
    ap.add_argument("--ntheta", type=int, default=72, help="number of angular bins")
    ap.add_argument("--times", default=None, metavar="T1,T2,...",
                    help="snapshot times (default: 3 times spread over the run)")
    ap.add_argument("--annulus", type=float, default=0.25,
                    help="snapshot annulus half-width relative to the front radius; "
                         "0 disables the restriction")
    ap.add_argument("--front-quantile", type=float, default=0.9,
                    help="energy quantile defining the front radius (as energy_rt.py)")
    ap.add_argument("--all-groups", action="store_true",
                    help="include the hybrid penalty groups 4,5 in the energy")
    ap.add_argument("--linear", action="store_true", help="linear instead of log color scale")
    args = ap.parse_args()

    t, centers, E_cell = read_cell_energies(args.file, args.all_groups)
    n_steps = len(t)

    if args.center is None:
        lo, hi = centers[:, :2].min(axis=0), centers[:, :2].max(axis=0)
        center = 0.5 * (lo + hi)
    else:
        center = np.array(list(map(float, args.center.split(","))))
    d = centers[:, :2] - center
    r = np.linalg.norm(d, axis=1)
    theta = np.arctan2(d[:, 1], d[:, 0])

    th_edges = np.linspace(-np.pi, np.pi, args.ntheta + 1)
    th_which = np.clip(np.digitize(theta, th_edges) - 1, 0, args.ntheta - 1)
    th_mid = 0.5 * (th_edges[:-1] + th_edges[1:])

    Eth = np.zeros((n_steps, args.ntheta))
    for s in range(n_steps):
        np.add.at(Eth[s], th_which, E_cell[s])

    # front radius per step (same definition as energy_rt.py, on unbinned radii)
    order = np.argsort(r)
    r_sorted = r[order]
    r_front = np.zeros(n_steps)
    for s in range(n_steps):
        cum = np.cumsum(E_cell[s][order])
        if cum[-1] > 0:
            r_front[s] = r_sorted[np.searchsorted(cum, args.front_quantile * cum[-1])]

    if args.times is None:
        idx = [max(1, n_steps // 4), n_steps // 2, (3 * n_steps) // 4]
    else:
        wanted = list(map(float, args.times.split(",")))
        idx = [int(np.argmin(np.abs(t - w))) for w in wanted]

    fig = plt.figure(figsize=(11, 4.2), constrained_layout=True)
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection="polar")

    pos = Eth[Eth > 0]
    norm = None if args.linear or pos.size == 0 else \
        LogNorm(vmin=max(pos.min(), pos.max() * 1e-6), vmax=pos.max())
    pc = ax1.pcolormesh(t, np.degrees(th_mid), Eth.T, norm=norm, cmap="inferno",
                        shading="nearest")
    fig.colorbar(pc, ax=ax1, label="energy per angular bin")
    ax1.set_xlabel("t"); ax1.set_ylabel(r"$\theta$ [deg]")
    ax1.set_yticks([-180, -90, 0, 90, 180])
    ax1.set_title(r"E($\theta$, t), all radii")

    # per-sector front radius: energy per angle is ~uniform even on a lattice, the
    # anisotropy lives in how FAR the front got per direction (grid: ~1.4x farther
    # along the axes than the diagonals, path-length zigzag factor sqrt(2))
    front_rows = []
    for s in idx:
        rf_theta = np.zeros(args.ntheta)
        for b in range(args.ntheta):
            m = th_which == b
            if not m.any() or E_cell[s][m].sum() <= 0:
                continue
            o = np.argsort(r[m])
            cum = np.cumsum(E_cell[s][m][o])
            rf_theta[b] = r[m][o][np.searchsorted(cum, args.front_quantile * cum[-1])]
        front_rows += [(t[s], np.degrees(th_mid[b]), rf_theta[b])
                       for b in range(args.ntheta)]
        ax2.plot(np.append(th_mid, th_mid[0]), np.append(rf_theta, rf_theta[0]),
                 lw=1, label=f"t={t[s]:.3g}")
    ax2.legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=8)
    ax2.set_title(f"front radius $r_{{{args.front_quantile:g}}}(\\theta)$")

    fig.savefig(args.save, dpi=150)
    print(f"wrote {args.save}")

    # pgfplots-ready sidecars (tikz versions read these)
    base = args.save.rsplit(".", 1)[0]
    with open(base + ".csv", "w") as f:
        f.write("t,theta_deg,E\n")
        for s in range(n_steps):
            for b in range(args.ntheta):
                f.write(f"{t[s]:.9g},{np.degrees(th_mid[b]):.6g},{Eth[s, b]:.9g}\n")
    with open(base + "-front.csv", "w") as f:
        f.write("t,theta_deg,r_front\n")
        for t_s, th_deg, rf in front_rows:
            f.write(f"{t_s:.9g},{th_deg:.6g},{rf:.9g}\n")
    print(f"wrote {base}.csv, {base}-front.csv")


if __name__ == "__main__":
    main()
