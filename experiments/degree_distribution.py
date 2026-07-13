#!/usr/bin/env python3
"""Compute the node degree distribution of a .geo.h5 hypergraph/network domain.

The domain is stored as a list of edges in /domain/edges of shape (n_edges, 2),
where each row holds the two endpoint node indices into /domain/points. The
degree of a node is the number of incident edges.

Usage:
    python3 degree_distribution.py output/test.geo.h5 [--plot] [--save out.png]
    python3 degree_distribution.py output/test.geo.h5 --spatial [--plot]
"""
import argparse

import h5py
import numpy as np


def load(path):
    """Return (degrees, points) where degrees[i] is the degree of node i
    and points[i] holds its (x, y, z) coordinates."""
    with h5py.File(path, "r") as f:
        edges = f["domain/edges"][:]          # (n_edges, 2) node indices
        points = f["domain/points"][:]        # (n_points, 3) coordinates
    # Each edge contributes +1 to the degree of both of its endpoints.
    degrees = np.bincount(edges.ravel(), minlength=points.shape[0])
    return degrees, points


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("geo_h5", help="path to a .geo.h5 file")
    p.add_argument("--spatial", action="store_true",
                   help="map mean degree over xy (hexbin) instead of a histogram")
    p.add_argument("--gridsize", type=int, default=80,
                   help="hexbin resolution for --spatial (default 80)")
    p.add_argument("--plot", action="store_true", help="show the plot via sxiv")
    p.add_argument("--save", metavar="PNG", help="save the plot to PNG")
    args = p.parse_args()

    degrees, points = load(args.geo_h5)
    n_points = points.shape[0]
    n_edges = degrees.sum() // 2

    print(f"file:      {args.geo_h5}")
    print(f"nodes:     {n_points}")
    print(f"edges:     {n_edges}")
    print(f"isolated:  {int(np.sum(degrees == 0))} nodes with degree 0")
    print(f"degree:    min={degrees.min()} max={degrees.max()} "
          f"mean={degrees.mean():.3f} median={int(np.median(degrees))}")
    print("\n degree   count    fraction")
    counts = np.bincount(degrees)
    for d, c in enumerate(counts):
        if c:
            print(f"   {d:4d}  {c:8d}   {c / n_points:8.4%}")

    if args.plot or args.save:
        import matplotlib
        if not args.plot:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if args.spatial:
            fig, ax = plt.subplots(figsize=(7, 6))
            x, y = points[:, 0], points[:, 1]
            hb = ax.hexbin(x, y, C=degrees, reduce_C_function=np.mean,
                           gridsize=args.gridsize, cmap="viridis")
            fig.colorbar(hb, ax=ax, label="mean degree")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_aspect("equal")
            ax.set_title(f"mean degree over xy: {args.geo_h5}")
            default_out = "/tmp/degree_spatial.png"
        else:
            fig, ax = plt.subplots()
            d = np.arange(len(counts))
            ax.bar(d[counts > 0], counts[counts > 0])
            ax.set_xlabel("degree")
            ax.set_ylabel("number of nodes")
            ax.set_yscale("log")
            ax.set_title(f"degree distribution: {args.geo_h5}")
            default_out = "/tmp/degree_distribution.png"
        fig.tight_layout()

        out = args.save or default_out
        fig.savefig(out, dpi=120)
        print(f"\nsaved plot to {out}")
        if args.plot:
            import subprocess
            subprocess.Popen(["sxiv", out])


if __name__ == "__main__":
    main()
