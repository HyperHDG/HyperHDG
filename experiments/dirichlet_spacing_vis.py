#!/usr/bin/env python
"""
dirichlet_spacing_vis.py -- visualize the Dirichlet-node spacing outliers that
set the boundary-density bound of R0 in gortz_constants.py.

For every Dirichlet node (types_points != 0) the spacing is its FULL-DIMENSION
distance to the nearest other Dirichlet node -- the same cKDTree k=2 query on
the same coordinates as gortz_constants.compute_R0, so the numbers here match
that report exactly.

The type code is decoded per make_geo.compute_types: base-3 with one digit per
axis (0 interior, 1 min face, 2 max face), type = d_x + 3 d_y + 9 d_z.  For a
thin 3D fiber slab the z faces are typically the sparsely-sampled ones, so the
outliers tend to live there.

Figure layout:
  left         : whole network in xy (edges, light gray), all Dirichlet nodes
                 colored by spacing (log scale), the --top worst circled with
                 their 3D gap radius drawn to scale;
  right top    : xy zoom around the outlier of --zoom-rank, with the nearest
                 Dirichlet neighbour and the empty gap disk;
  right bottom : xz view of the same window (full z range), showing where in
                 the slab thickness the outlier and its neighbours sit.

A table of the --top worst nodes (index, xyz, decoded faces, spacing) is
printed.

Usage:
    python dirichlet_spacing_vis.py domain.geo.h5
    python dirichlet_spacing_vis.py domain.geo.h5 -o spacing.png --top 10
    python dirichlet_spacing_vis.py domain.geo.h5 --zoom-rank 2
"""

import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree


def load_domain(path):
    with h5py.File(path, "r") as f:
        g = f["domain"]
        points = g["points"][:].astype(float)
        edges = g["edges"][:].astype(np.int64)
        types_points = None
        if "types_points" in g:
            types_points = np.asarray(g["types_points"][:]).reshape(-1)
    return points, edges, types_points


def face_label(t):
    """Decode a make_geo base-3 type code into face names like 'z-' or 'x+y-'."""
    parts = []
    for axis in "xyz":
        d = t % 3
        t //= 3
        if d == 1:
            parts.append(axis + "-")
        elif d == 2:
            parts.append(axis + "+")
    return "".join(parts) or "int"


def parse_types(spec):
    """'1-10' or '9,18' -> set of type codes; None = all nonzero (gortz proxy)."""
    if spec is None:
        return None
    out = set()
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def dirichlet_spacing(points, types_points, types=None):
    """Nearest-neighbour spacing among Dirichlet nodes, exactly as compute_R0:
    full-dimensional coordinates, k=2 query (k=1 is self).  `types` restricts
    which type codes count as Dirichlet (e.g. the solver's clamp set 1..10);
    None matches gortz_constants (all nonzero)."""
    if types is None:
        mask = types_points != 0
    else:
        mask = np.isin(types_points, list(types))
    didx = np.nonzero(mask)[0]
    dpts = points[didx]
    tree = cKDTree(dpts)
    d, nn = tree.query(dpts, k=2)
    return didx, dpts, d[:, 1], nn[:, 1]


def draw_edges(ax, points, edges, cols=(0, 1), lw=0.15, color="0.75"):
    segs = points[:, cols][edges]    # (n_edges, 2, 2)
    ax.add_collection(LineCollection(segs, linewidths=lw, colors=color,
                                     rasterized=True, zorder=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("domain", help="input .geo.h5 domain (make_geo2 format)")
    ap.add_argument("-o", "--output", default="dirichlet-spacing.png",
                    help="output image (default dirichlet-spacing.png)")
    ap.add_argument("--top", type=int, default=10,
                    help="number of worst-spacing nodes to mark/list (default 10)")
    ap.add_argument("--zoom-rank", type=int, default=1,
                    help="which outlier the zoom panels center on, 1 = worst")
    ap.add_argument("--zoom-halfwidth", type=float, default=None,
                    help="half-width of the zoom window (default 2x the gap)")
    ap.add_argument("--types", default=None,
                    help="type codes counted as Dirichlet, e.g. '1-10' for the "
                         "solver clamp set (default: all nonzero, as in "
                         "gortz_constants)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    points, edges, types_points = load_domain(args.domain)
    if types_points is None:
        raise SystemExit("file has no types_points dataset -- no Dirichlet nodes")

    didx, dpts, spacing, nn_local = dirichlet_spacing(points, types_points,
                                                      parse_types(args.types))
    order = np.argsort(spacing)[::-1]
    top = order[:args.top]

    print(f"{len(didx)} Dirichlet nodes; spacing p50={np.percentile(spacing, 50):.4g} "
          f"p99={np.percentile(spacing, 99):.4g} max={spacing.max():.4g}")
    print(f"{'rank':>4} {'node':>8} {'x':>10} {'y':>10} {'z':>8} "
          f"{'faces':>6} {'spacing':>9}")
    for r, k in enumerate(top, 1):
        print(f"{r:>4} {didx[k]:>8} {dpts[k, 0]:>10.2f} {dpts[k, 1]:>10.2f} "
              f"{dpts[k, 2]:>8.2f} {face_label(int(types_points[didx[k]])):>6} "
              f"{spacing[k]:>9.4g}")

    fig = plt.figure(figsize=(17, 8.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[2.2, 1])
    ax0 = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 1])

    # ---- left: whole domain, xy projection --------------------------------
    draw_edges(ax0, points, edges)
    o = np.argsort(spacing)          # draw big-spacing nodes last (on top)
    sc = ax0.scatter(dpts[o, 0], dpts[o, 1], c=spacing[o], s=2,
                     norm=LogNorm(vmin=max(spacing.min(), 1e-3),
                                  vmax=spacing.max()),
                     cmap="viridis", zorder=2)
    fig.colorbar(sc, ax=ax0, label="nearest-Dirichlet-neighbour spacing (3D)",
                 shrink=0.85)
    for r, k in enumerate(top, 1):
        ax0.add_patch(plt.Circle(dpts[k, :2], spacing[k], fill=False,
                                 color="red", lw=1.0, zorder=3))
        ax0.annotate(f"#{r}: {spacing[k]:.3g} ({face_label(int(types_points[didx[k]]))})",
                     dpts[k, :2], xytext=(6, 6), textcoords="offset points",
                     color="red", fontsize=8, zorder=4)
    ax0.set_title(f"Dirichlet-node spacing, top {args.top} circled (3D gap to scale)")
    ax0.set_aspect("equal")
    ax0.autoscale_view()

    # ---- right: zoom on one outlier, xy and xz -----------------------------
    k = order[args.zoom_rank - 1]
    c = dpts[k]
    nb = dpts[nn_local[k]]
    gap = spacing[k]
    hw = args.zoom_halfwidth if args.zoom_halfwidth is not None else 2.0 * gap
    in_win = ((np.abs(dpts[:, 0] - c[0]) <= hw) &
              (np.abs(dpts[:, 1] - c[1]) <= hw))
    # edges with at least one endpoint in the xy window, so the xz projection
    # shows only the local fibers instead of the whole domain collapsed in y
    node_in = ((np.abs(points[:, 0] - c[0]) <= hw) &
               (np.abs(points[:, 1] - c[1]) <= hw))
    eloc = edges[node_in[edges[:, 0]] | node_in[edges[:, 1]]]

    for ax, cols in ((ax1, (0, 1)), (ax2, (0, 2))):
        draw_edges(ax, points, eloc, cols=cols, lw=0.4, color="0.6")
        ax.scatter(dpts[in_win, cols[0]], dpts[in_win, cols[1]], s=20,
                   color="tab:blue", zorder=3, label="Dirichlet nodes")
        ax.scatter(c[cols[0]], c[cols[1]], s=110, color="red", zorder=4,
                   marker="*",
                   label=f"outlier #{args.zoom_rank} (node {didx[k]})")
        ax.scatter(nb[cols[0]], nb[cols[1]], s=45, color="orange", zorder=4,
                   label="nearest Dirichlet neighbour")
        ax.plot([c[cols[0]], nb[cols[0]]], [c[cols[1]], nb[cols[1]]],
                "r--", lw=1, zorder=3)
        ax.set_xlim(c[0] - hw, c[0] + hw)

    ax1.add_patch(plt.Circle(c[:2], gap, fill=False, color="red", lw=1.2,
                             ls="--", zorder=3))
    ax1.set_ylim(c[1] - hw, c[1] + hw)
    ax1.set_aspect("equal")
    ax1.set_title(f"xy zoom: node {didx[k]} at ({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f}), "
                  f"faces {face_label(int(types_points[didx[k]]))}, gap = {gap:.4g}")
    ax1.legend(loc="upper right", fontsize=8)

    zlo, zhi = points[:, 2].min(), points[:, 2].max()
    pad = 0.05 * (zhi - zlo)
    ax2.set_ylim(zlo - pad, zhi + pad)
    ax2.set_title("same window, xz view (full slab thickness)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")

    fig.suptitle(args.domain)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
