#!/usr/bin/env python3
"""Read per-cell, per-step energies from a wave.vtkhdf and plot totals over time.

Components in CellData/energies are laid out (block of size space_dim each):
  0: ½ ∫ n²/Cn        (axial strain)
  1: ½ ∫ m²/Cm        (bending strain)
  2: ½ ∫ v²/Cu        (translational kinetic)
  3: ½ ∫ s²/Cr        (rotational kinetic)
  4: ½ τ Σ_bdr ∫ (u-λu)²   (hybrid translational)
  5: ½ τ Σ_bdr ∫ (r-λr)²   (hybrid rotational)

Within each block the d-th entry is the contribution from edge-local direction d.
We sum over directions to get 6 physical energies, and over cells to get a total per step.
"""

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


GROUP_LABELS = [
    r"$\frac{1}{2}\int n^2/C_n$  axial strain",
    r"$\frac{1}{2}\int m^2/C_m$  bending strain",
    r"$\frac{1}{2}\int v^2/C_u$  trans. kinetic",
    r"$\frac{1}{2}\int s^2/C_r$  rot. kinetic",
    r"$\frac{1}{2}\tau\,\sum\!\int(u-\lambda_u)^2$  hybrid trans.",
    r"$\frac{1}{2}\tau\,\sum\!\int(r-\lambda_r)^2$  hybrid rot.",
]


def read_energies(path):
    with h5py.File(path, "r") as f:
        E = f["/VTKHDF/CellData/energies"][...]              # (n_steps*n_cells, n_comps)
        t = f["/VTKHDF/Steps/Values"][...]                   # (n_steps,)
        off = f["/VTKHDF/Steps/CellDataOffsets/energies"][...]  # (n_steps,) row offsets
    n_steps = t.shape[0]
    n_comps = E.shape[1]
    assert n_comps % 6 == 0, f"expected n_comps divisible by 6, got {n_comps}"
    space_dim = n_comps // 6

    # rows for step s: [off[s], off[s+1]); last step runs to end
    bounds = np.append(off, E.shape[0])
    per_step = np.zeros((n_steps, n_comps), dtype=E.dtype)
    for s in range(n_steps):
        per_step[s] = E[bounds[s]:bounds[s+1]].sum(axis=0)

    # collapse the d-axis: (n_steps, 6, space_dim) -> (n_steps, 6)
    grouped = per_step.reshape(n_steps, 6, space_dim).sum(axis=2)
    return t, grouped, per_step, space_dim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?", default="output/energy/wave.vtkhdf",
                    help="path to wave.vtkhdf")
    ap.add_argument("-o", "--save", help="save plot to this path instead of showing")
    ap.add_argument("--log", action="store_true", help="log scale on y")
    ap.add_argument("--all", action="store_true",
                    help="plot all 6*space_dim raw components instead of grouped")
    args = ap.parse_args()

    t, grouped, raw, space_dim = read_energies(args.file)
    total = grouped.sum(axis=1)

    physical = grouped[:, :4].sum(axis=1)  # n, m strain + v, s kinetic
    hybrid   = grouped[:, 4:].sum(axis=1)  # u, r hybrid penalty (HDG stabilization)

    fig, ax = plt.subplots(figsize=(8, 5))
    if args.all:
        for c in range(raw.shape[1]):
            g, d = divmod(c, space_dim)
            ax.plot(t, raw[:, c], label=f"{GROUP_LABELS[g].split('  ')[1]} (d={d})",
                    linewidth=1)
    else:
        for c in range(grouped.shape[1]):
            ax.plot(t, grouped[:, c], label=GROUP_LABELS[c], linewidth=1)
    ax.plot(t, physical, label="physical (n+m+v+s)", color="tab:green", linewidth=2, linestyle="--")
    ax.plot(t, hybrid,   label="hybrid (u+r, HDG stab.)", color="tab:red", linewidth=2, linestyle="--")
    ax.plot(t, total, label="total", color="k", linewidth=2)

    ax.set_xlabel("time")
    ax.set_ylabel("energy")
    if args.log:
        ax.set_yscale("log")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if args.save:
        Path(args.save).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.save, bbox_inches="tight", pad_inches=0.05)
        print(f"wrote {args.save}", file=sys.stderr)
    else:
        plt.show()


if __name__ == "__main__":
    main()
