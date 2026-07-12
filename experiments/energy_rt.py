#!/usr/bin/env python3
"""Radial energy transport E(r,t) from a wave.vtkhdf with per-cell energies.

Reads CellData/energies (layout as in energy.py: 6 groups x space_dim, see there),
bins the per-cell energy over the distance r of the (static) cell center from the tap
center, and plots
  (1) the E(r,t) heatmap -- a ballistic front is a straight line of slope c, scattering
      broadens it, diffusive transport bends it,
  (2) the energy-weighted mean square radius <r^2>(t) with t and t^2 guide lines
      (ballistic <r^2> ~ t^2, diffusive ~ t), and
  (3) the front radius r_q(t) (radius containing --front-quantile of the energy).

By default only the physical energy groups 0..3 (axial/bending strain, translational/
rotational kinetic) are summed; --all-groups adds the hybrid penalty terms.

This is the post-hoc fallback for small/medium runs; the in-situ histograms written by
the solver (upcoming) produce the same binning without storing per-cell energies.
"""

import argparse

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def read_cell_energies(path, all_groups=False):
    """Return (t, centers, E_cell) with E_cell of shape (n_steps, n_cells)."""
    with h5py.File(path, "r") as f:
        r = f["VTKHDF"]
        E = r["CellData/energies"][...]
        t = r["Steps/Values"][...]
        off = r["Steps/CellDataOffsets/energies"][...]
        pts = r["Points"][...]
        conn = r["Connectivity"][...]
    n_steps = t.shape[0]
    n_comps = E.shape[1]
    assert n_comps % 6 == 0, f"expected components divisible by 6, got {n_comps}"
    # files written before the timowave energy() guard carry NaN kinetic entries on
    # massless welds (0/0); physically their energy is zero
    n_bad = np.isnan(E).sum() + np.isinf(E).sum()
    if n_bad:
        print(f"note: {n_bad} NaN/Inf energy entries (massless welds) treated as 0")
        E = np.nan_to_num(E, nan=0.0, posinf=0.0, neginf=0.0)
    space_dim = n_comps // 6
    n_groups = 6 if all_groups else 4
    bounds = np.append(off, E.shape[0])
    n_cells = int(bounds[1] - bounds[0])

    # physical energy per cell: sum selected groups over local directions
    sel = E[:, : n_groups * space_dim].sum(axis=1)
    E_cell = np.zeros((n_steps, n_cells))
    for s in range(n_steps):
        E_cell[s] = sel[bounds[s]:bounds[s + 1]]

    centers = 0.5 * (pts[conn[0::2]] + pts[conn[1::2]])  # static line cells
    assert centers.shape[0] == n_cells
    return t, centers, E_cell


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file", help="path to wave.vtkhdf (needs CellData/energies)")
    ap.add_argument("-o", "--save", default="energy_rt.png", help="output image path")
    ap.add_argument("--center", default=None, metavar="X,Y",
                    help="tap center (default: xy bbox center of the mesh)")
    ap.add_argument("--nr", type=int, default=128, help="number of radial bins")
    ap.add_argument("--all-groups", action="store_true",
                    help="include the hybrid penalty groups 4,5 in the energy")
    ap.add_argument("--front-quantile", type=float, default=0.9,
                    help="energy quantile defining the front radius")
    ap.add_argument("--linear", action="store_true", help="linear instead of log color scale")
    args = ap.parse_args()

    t, centers, E_cell = read_cell_energies(args.file, args.all_groups)

    if args.center is None:
        lo, hi = centers[:, :2].min(axis=0), centers[:, :2].max(axis=0)
        center = 0.5 * (lo + hi)
    else:
        center = np.array(list(map(float, args.center.split(","))))
    r = np.linalg.norm(centers[:, :2] - center, axis=1)

    r_edges = np.linspace(0.0, r.max(), args.nr + 1)
    which = np.clip(np.digitize(r, r_edges) - 1, 0, args.nr - 1)

    n_steps = len(t)
    Ert = np.zeros((n_steps, args.nr))
    for s in range(n_steps):
        np.add.at(Ert[s], which, E_cell[s])
    # bins containing no cells are undefined, not zero-energy -- mask them out
    count = np.bincount(which, minlength=args.nr)
    Ert[:, count == 0] = np.nan

    tot = np.nansum(Ert, axis=1)
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        msr = np.nansum(Ert * r_mid**2, axis=1) / tot  # <r^2>(t)
        cum = np.nancumsum(Ert, axis=1) / np.maximum(tot, 1e-300)[:, None]
    r_front = np.array([np.interp(args.front_quantile, cum[s], r_mid) for s in range(n_steps)])

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11, 4.2), width_ratios=[1.6, 1], constrained_layout=True)

    with np.errstate(invalid="ignore"):
        pos = Ert[np.nan_to_num(Ert) > 0]
    norm = None if args.linear or pos.size == 0 else \
        LogNorm(vmin=max(pos.min(), pos.max() * 1e-8), vmax=pos.max())
    pc = ax1.pcolormesh(t, r_mid, Ert.T, norm=norm, cmap="inferno", shading="nearest")
    ax1.plot(t, r_front, color="cyan", lw=1,
             label=f"$r_{{{args.front_quantile:g}}}(t)$ front")
    fig.colorbar(pc, ax=ax1, label="energy per radial bin")
    ax1.set_xlabel("t"); ax1.set_ylabel("r"); ax1.legend(loc="upper left")
    ax1.set_title("E(r, t)")

    m = tot > 0
    ax2.loglog(t[m], msr[m], label=r"$\langle r^2\rangle(t)$")
    if m.sum() >= 3:
        t0 = t[m][min(max(1, len(t[m]) // 10), len(t[m]) - 1)]
        y0 = np.interp(t0, t[m], msr[m])
        ax2.loglog(t[m], y0 * (t[m] / t0) ** 2, "k--", lw=0.8,
                   label=r"$\propto t^2$ (ballistic)")
        ax2.loglog(t[m], y0 * (t[m] / t0), "k:", lw=0.8, label=r"$\propto t$ (diffusive)")
    ax2.set_xlabel("t"); ax2.set_ylabel(r"$\langle r^2\rangle$"); ax2.legend()
    ax2.set_title("energy spread")

    fig.savefig(args.save, dpi=150)
    print(f"wrote {args.save}")

    # pgfplots-ready sidecars (tikz versions read these): long-format E(r,t) for a
    # matrix/surf plot, and the spread/front curves
    base = args.save.rsplit(".", 1)[0]
    with open(base + ".csv", "w") as f:
        f.write("t,r,E\n")
        for s in range(n_steps):
            for j in range(args.nr):
                if not np.isnan(Ert[s, j]):
                    f.write(f"{t[s]:.9g},{r_mid[j]:.9g},{Ert[s, j]:.9g}\n")
    with open(base + "-spread.csv", "w") as f:
        f.write("t,total,msr,r_front\n")
        for s in range(n_steps):
            f.write(f"{t[s]:.9g},{tot[s]:.9g},{msr[s]:.9g},{r_front[s]:.9g}\n")
    print(f"wrote {base}.csv, {base}-spread.csv")


if __name__ == "__main__":
    main()
