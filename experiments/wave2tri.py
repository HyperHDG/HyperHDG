#!/usr/bin/env python3
"""Convert a broken-space network wave.vtkhdf into a conforming triangulated one.

The plotter writes the non-conforming (per-edge) solution: every edge carries its own
copies of its endpoints, so a node appears deg(v) times and PointData/values holds the
edge-local solutions (discontinuous across nodes at HDG-jump level). For planar grids
that makes surface rendering absurdly expensive: netvis --surface must merge points and
Delaunay-triangulate INSIDE the render loop, re-executed for every frame (and twice per
frame with --color-rescale frame).

This converter does the geometry work exactly once, offline:
  1. merge bit-exact duplicate points (grid endpoint coords are identical by
     construction) -> the true nodes,
  2. Delaunay-triangulate the nodes' xy once (--alpha bounds the circumradius so
     concave domains/holes are not bridged; 0 = convex hull, fine for grids),
  3. average the per-edge values onto nodes per frame (the duplicates differ only by
     the inter-edge jump, so the average is the natural conforming/trace-like field),
  4. write a static-mesh transient VTKHDF of triangles.

The output renders with plain netvis (no --surface): warp + color on a triangle sheet.
"""

import argparse

import h5py
import numpy as np
from scipy.spatial import Delaunay


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input", help="broken-space wave.vtkhdf (line cells)")
    ap.add_argument("output", help="triangulated conforming wave.vtkhdf to write")
    ap.add_argument("--alpha", type=float, default=0.,
                    help="drop triangles with circumradius above ALPHA (~2x mean edge "
                         "length keeps holes open in disordered nets); 0 = convex hull")
    args = ap.parse_args()

    with h5py.File(args.input, "r") as f:
        r = f["VTKHDF"]
        pts = r["Points"][...]
        vals = r["PointData/values"]
        t = r["Steps/Values"][...]
        v_off = r["Steps/PointDataOffsets/values"][...]
        n_comp = vals.shape[1]
        n_pts = pts.shape[0]
        n_steps = len(t)

        # 1. merge bit-exact duplicates
        _, uniq_idx, inv = np.unique(pts.view([("", pts.dtype)] * 3), axis=0,
                                     return_index=True, return_inverse=True)
        nodes = pts[uniq_idx]
        count = np.bincount(inv).astype(vals.dtype)
        n_nodes = nodes.shape[0]
        print(f"{n_pts} plot points -> {n_nodes} nodes "
              f"(mean multiplicity {n_pts / n_nodes:.2f})")

        # 2. triangulate once
        tri = Delaunay(nodes[:, :2].astype(np.float64))
        cells = tri.simplices
        if args.alpha > 0:
            a, b, c = (tri.points[cells[:, i]] for i in range(3))
            la, lb, lc = (np.linalg.norm(x - y, axis=1) for x, y in ((b, c), (a, c), (a, b)))
            area = 0.5 * np.abs(np.cross(b - a, c - a))
            circum = la * lb * lc / np.maximum(4 * area, 1e-300)
            cells = cells[circum < args.alpha]
        n_cells = cells.shape[0]
        print(f"{n_cells} triangles")

        # 4. output skeleton (static mesh, extendable values)
        with h5py.File(args.output, "w") as g:
            root = g.create_group("VTKHDF")
            root.attrs.create("Version", [2, 0], dtype="int64")
            root.attrs.create("Type", np.bytes_("UnstructuredGrid"))
            root.create_dataset("Points", data=nodes)
            root.create_dataset("Connectivity", data=cells.ravel().astype(np.int64))
            root.create_dataset("Offsets",
                                data=np.arange(0, 3 * n_cells + 3, 3, dtype=np.int64))
            root.create_dataset("Types", data=np.full(n_cells, 5, dtype=np.uint8))
            root.create_dataset("NumberOfPoints", data=np.array([n_nodes], dtype=np.int64))
            root.create_dataset("NumberOfCells", data=np.array([n_cells], dtype=np.int64))
            root.create_dataset("NumberOfConnectivityIds",
                                data=np.array([3 * n_cells], dtype=np.int64))
            pd = root.create_group("PointData")
            out_vals = pd.create_dataset("values", shape=(n_steps * n_nodes, n_comp),
                                         dtype=vals.dtype)

            # 3. per-frame fold: average edge-local values onto nodes
            for s in range(n_steps):
                block = vals[v_off[s]:v_off[s] + n_pts]
                acc = np.zeros((n_nodes, n_comp), dtype=np.float64)
                np.add.at(acc, inv, block.astype(np.float64))
                out_vals[s * n_nodes:(s + 1) * n_nodes] = (acc / count[:, None])

            steps = root.create_group("Steps")
            steps.attrs.create("NSteps", n_steps, dtype="int64")
            steps.create_dataset("Values", data=t)
            zeros = np.zeros(n_steps, dtype=np.int64)
            steps.create_dataset("PartOffsets", data=zeros)
            steps.create_dataset("PointOffsets", data=zeros)
            steps.create_dataset("CellOffsets", data=zeros)
            steps.create_dataset("ConnectivityIdOffsets", data=zeros)
            steps.create_dataset("NumberOfParts", data=np.ones(n_steps, dtype=np.int64))
            pdo = steps.create_group("PointDataOffsets")
            pdo.create_dataset("values",
                               data=np.arange(n_steps, dtype=np.int64) * n_nodes)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
