#!/usr/bin/env python
"""
gortz_constants.py -- estimate the network constants of Assumption 3.5 in

  M. Goertz, F. Hellman, A. Maalqvist,
  "Iterative solution of spatial network models by subspace decomposition",
  Math. Comp. 93 (2024), pp. 231-258.  (gortz.pdf, Section 6.1, pp. 250-253)

for a domain in .geo.h5 format as produced by experiments/make_geo2.py.

Two constants are reported (Assumption 3.5, p. 242):

  sigma  (homogeneity, 3.5.1)
      max_{B_R(x) c Omega} |1|^2_{M,B_R(x)}  <=  sigma * min_{...} |1|^2_{M,B_R(x)}
      estimated on p. 251 by
          sigma(R) = max_{T in T_H} |1|^2_{M,T}  /  min_{T in T_H} |1|^2_{M,T}
      over the cells T of a regular n x n grid of side R = 1/n on the domain.
      Here the M-norm of the constant function 1 restricted to a cell is the
      network mass in that cell:
          |1|^2_{M,T} = sum_{x in T} m_x ,   m_x = (1/2) sum_{e ~ x} mass_e
      i.e. every edge contributes half its mass to each of its two endpoints
      (the row sum of the consistent/lumped mass matrix M).  For density 1 the
      edge mass equals the edge length, which is the paper's convention
      (Section 6.1: "mass is defined as the total edge length").

  R0  (length scale)
      The smallest scale at which the four network assumptions can hold.  Two
      of them give a directly computable geometric lower bound:
          locality        (3.5.3):  R0 > max_e |e|          (longest edge)
          boundary density (3.5.4):  R0 > max_{y in Gamma} dist(y, N(Gamma))
      We report the longest edge (the dominant, unambiguous scale) and, as a
      proxy for boundary density, the largest nearest-neighbour spacing between
      Dirichlet nodes (types_points != 0).  R0 = max of the two.

Optionally (--mu) the connectivity constant mu of Lemma 3.6 is estimated per
cell from the generalized eigenproblem  Lbar u = lambda_2 Mbar u  (Neumann),
mu(x,R) = R^-1 lambda_2^-1/2,  mu = max_{x} mu(x,R).  This is the more expensive
and more approximate quantity; see notes in estimate_mu().

Usage:
    python gortz_constants.py domain.geo.h5
    python gortz_constants.py domain.geo.h5 --cells 4 8 16 32 64
    python gortz_constants.py domain.geo.h5 --use-properties-mass
    python gortz_constants.py domain.geo.h5 --mu --cells 8 16 32
"""

import argparse
import time

import h5py
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree


def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)


def load_domain(path):
    """Read points, edges, optional per-edge mass (properties col 0) and the
    Dirichlet marker types_points from a make_geo2 .geo.h5 file."""
    with h5py.File(path, "r") as f:
        g = f["domain"]
        points = g["points"][:].astype(float)
        edges = g["edges"][:].astype(np.int64)
        mass = g["properties"][:, 0].astype(float) if "properties" in g else None
        types_points = None
        if "types_points" in g:
            types_points = np.asarray(g["types_points"][:]).reshape(-1)
        size = np.array(g.attrs["size"]) if "size" in g.attrs else None
    tprint(f"loaded '{path}': {len(points)} nodes, {len(edges)} edges, "
           f"properties={'yes' if mass is not None else 'no'}")
    return dict(points=points, edges=edges, mass=mass,
                types_points=types_points, size=size)


def edge_lengths(points, edges):
    return np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)


def edge_mass(dom, use_properties_mass, lengths):
    """Per-edge mass used for the M-norm.  Default is the geometric edge length
    (density 1, the paper's convention and always strictly positive); with
    --use-properties-mass we take properties column 0, the mass the solver's M
    actually uses (may contain zeros for virtual/zero-mass fibers)."""
    if use_properties_mass:
        if dom["mass"] is None:
            raise SystemExit("--use-properties-mass: file has no 'properties' dataset")
        m = dom["mass"].copy()
        n_zero = int((m <= 0).sum())
        if n_zero:
            tprint(f"  note: {n_zero}/{len(m)} edges have non-positive properties mass "
                   f"(virtual/zero-mass fibers); cells containing only such edges "
                   f"count as empty for sigma")
        return m
    return lengths


def compute_R0(points, edges, types_points, lengths):
    """Geometric length scale R0 (Assumption 3.5.3-4).

    Returns (R0, info dict).  R0 = max(longest edge, boundary-density proxy)."""
    l_max = float(lengths.max())
    info = {"l_max": l_max,
            "l_pctl": {p: float(np.percentile(lengths, p))
                       for p in (50, 90, 99, 100)}}

    bnd = None
    if types_points is not None:
        dnodes = points[types_points != 0]
        if len(dnodes) > 1:
            # covering proxy: nearest-neighbour spacing among Dirichlet nodes.
            # The continuous Dirichlet boundary Gamma is covered by these nodes,
            # so the largest gap between neighbours bounds dist(y, N(Gamma)).
            tree = cKDTree(dnodes)
            d, _ = tree.query(dnodes, k=2)   # k=1 is self (distance 0)
            spacing = d[:, 1]
            bnd = float(spacing.max())
            info["n_dirichlet"] = int(len(dnodes))
            info["bnd_spacing_max"] = bnd
            info["bnd_spacing_pctl"] = {p: float(np.percentile(spacing, p))
                                        for p in (50, 90, 99, 100)}

    R0 = l_max if bnd is None else max(l_max, bnd)
    info["R0"] = R0
    info["R0_from"] = "longest-edge" if (bnd is None or l_max >= bnd) else "boundary-density"
    return R0, info


def sigma_grid(points, edges, m_edge, n):
    """sigma(R) on a regular n x n tiling (side R = extent/n) of the xy bbox.

    |1|^2_{M,T} = sum over nodes in cell T of the node mass, where each edge
    donates half its mass to each endpoint.  sigma = max cell / min positive
    cell.  Returns a stats dict."""
    xy = points[:, :2]
    lo = xy.min(0)
    ext = xy.max(0) - lo
    ext = np.where(ext > 0, ext, 1.0)          # guard degenerate (1-D) axis

    ix = np.clip(((xy[:, 0] - lo[0]) / ext[0] * n).astype(int), 0, n - 1)
    iy = np.clip(((xy[:, 1] - lo[1]) / ext[1] * n).astype(int), 0, n - 1)
    cell = iy * n + ix

    node_mass = np.zeros(len(points))
    np.add.at(node_mass, edges[:, 0], m_edge / 2.0)
    np.add.at(node_mass, edges[:, 1], m_edge / 2.0)

    cell_mass = np.bincount(cell, weights=node_mass, minlength=n * n)
    pos = cell_mass[cell_mass > 0]
    n_empty = int((cell_mass == 0).sum())
    sigma = float(cell_mass.max() / pos.min()) if pos.size else float("inf")
    return dict(n=n, R=float(ext.max() / n), sigma=sigma,
                cmin=float(pos.min()) if pos.size else 0.0,
                cmax=float(cell_mass.max()),
                mean=float(cell_mass.mean()), n_cells=n * n, n_empty=n_empty)


def estimate_mu(points, edges, m_edge, lengths, R0, n, max_nodes=40000):
    """Connectivity constant mu(R) of Lemma 3.6, Neumann case, per cell.

    For each grid cell B_R(x) we form the induced subgraph on the nodes inside
    the R0-enlarged cell B_{R+R0}(x), build the weighted graph Laplacian Lbar
    (edge weight 1/|e|, matching operator L of Example 3.2 with gamma=1) and the
    lumped mass Mbar, and take lambda_2 = smallest nonzero eigenvalue of
    Lbar u = lambda Mbar u.  Then mu(x,R) = R^-1 lambda_2^-1/2 and mu = max.

    This is an ESTIMATE: it uses the induced subgraph on the enlarged cell
    rather than the BFS-minimal connected subgraph of the paper, and cells whose
    subgraph exceeds max_nodes are skipped (reported).  Disconnected subgraphs
    (lambda_2 ~ 0) signal a connectivity violation and are reported, not
    silently dropped."""
    from scipy.sparse.linalg import eigsh

    xy = points[:, :2]
    lo = xy.min(0)
    ext = xy.max(0) - lo
    ext = np.where(ext > 0, ext, 1.0)
    hx, hy = ext / n
    R = float(ext.max() / n)

    w = 1.0 / lengths                              # Laplacian edge weights
    node_mass = np.zeros(len(points))
    np.add.at(node_mass, edges[:, 0], m_edge / 2.0)
    np.add.at(node_mass, edges[:, 1], m_edge / 2.0)

    mu_max, worst_cell = 0.0, None
    n_skip, n_disc, n_eval = 0, 0, 0
    for i in range(n):
        for j in range(n):
            x0, x1 = lo[0] + i * hx, lo[0] + (i + 1) * hx
            y0, y1 = lo[1] + j * hy, lo[1] + (j + 1) * hy
            # nodes of the R0-enlarged cell
            sel = ((xy[:, 0] >= x0 - R0) & (xy[:, 0] <= x1 + R0) &
                   (xy[:, 1] >= y0 - R0) & (xy[:, 1] <= y1 + R0))
            idx = np.nonzero(sel)[0]
            if idx.size < 2:
                continue
            if idx.size > max_nodes:
                n_skip += 1
                continue

            remap = -np.ones(len(points), dtype=np.int64)
            remap[idx] = np.arange(idx.size)
            emask = sel[edges[:, 0]] & sel[edges[:, 1]]
            se = remap[edges[emask]]
            if se.size == 0:
                continue
            we = w[emask]

            # weighted Laplacian L = D - W  (symmetric)
            a, b = se[:, 0], se[:, 1]
            W = sp.coo_matrix((we, (a, b)), shape=(idx.size, idx.size))
            W = (W + W.T)
            L = sp.diags(np.asarray(W.sum(1)).ravel()) - W
            M = sp.diags(node_mass[idx])
            if node_mass[idx].min() <= 0:
                continue                            # ill-defined M

            try:
                # smallest two eigenvalues of the SPD pencil; shift just below 0
                vals = eigsh(L.tocsc(), k=2, M=M.tocsc(),
                             sigma=-1e-12 * (we.sum()), which="LM",
                             return_eigenvectors=False)
                vals = np.sort(vals)
                lam2 = vals[1]
            except Exception:
                continue

            n_eval += 1
            scale = we.mean()
            if lam2 <= 1e-10 * scale:               # subgraph disconnected
                n_disc += 1
                continue
            mu_cell = (1.0 / R) * lam2 ** -0.5
            if mu_cell > mu_max:
                mu_max, worst_cell = mu_cell, (i, j)

    return dict(n=n, R=R, mu=mu_max, worst_cell=worst_cell,
                n_eval=n_eval, n_skip=n_skip, n_disconnected=n_disc)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("domain", help="input .geo.h5 domain (make_geo2 format)")
    ap.add_argument("--cells", type=int, nargs="+", default=[4, 8, 16, 32, 64],
                    metavar="N",
                    help="cells per axis for the sigma sweep (R = extent/N); "
                         "default reproduces the paper's R^-1 = 4..64")
    ap.add_argument("--use-properties-mass", action="store_true",
                    help="use properties column 0 as edge mass instead of the "
                         "geometric edge length (default)")
    ap.add_argument("--mu", action="store_true",
                    help="also estimate the connectivity constant mu (slow)")
    ap.add_argument("--mu-max-nodes", type=int, default=40000,
                    help="skip cells whose enlarged subgraph exceeds this many "
                         "nodes in the mu estimate (default 40000)")
    args = ap.parse_args()

    dom = load_domain(args.domain)
    lengths = edge_lengths(dom["points"], dom["edges"])
    m_edge = edge_mass(dom, args.use_properties_mass, lengths)

    xy_ext = (dom["points"][:, :2].max(0) - dom["points"][:, :2].min(0))
    tprint(f"xy bounding box extent: {xy_ext[0]:.4g} x {xy_ext[1]:.4g}")

    # ---- R0 --------------------------------------------------------------
    R0, r0i = compute_R0(dom["points"], dom["edges"], dom["types_points"], lengths)
    print()
    print("=== R0  (length scale, Assumption 3.5.3-4) ===")
    lp = r0i["l_pctl"]
    print(f"  longest edge (locality bound)        : {r0i['l_max']:.4g}")
    print(f"    edge-length pctl  p50={lp[50]:.4g}  p90={lp[90]:.4g}  "
          f"p99={lp[99]:.4g}  max={lp[100]:.4g}")
    if "bnd_spacing_max" in r0i:
        bp = r0i["bnd_spacing_pctl"]
        print(f"  Dirichlet-node spacing (bnd bound)   : {r0i['bnd_spacing_max']:.4g}"
              f"  ({r0i['n_dirichlet']} Dirichlet nodes)")
        print(f"    spacing pctl      p50={bp[50]:.4g}  p90={bp[90]:.4g}  "
              f"p99={bp[99]:.4g}  max={bp[100]:.4g}")
    else:
        print("  Dirichlet-node spacing               : (no types_points / <2 nodes)")
    print(f"  --> R0 = {R0:.4g}   (set by {r0i['R0_from']})")
    print(f"      relative to xy extent: R0 / max_extent = {R0 / xy_ext.max():.4g}"
          f"   (i.e. R0^-1 ~ {xy_ext.max() / R0:.1f})")

    # ---- sigma -----------------------------------------------------------
    massname = "properties col 0" if args.use_properties_mass else "edge length"
    print()
    print(f"=== sigma  (homogeneity, Assumption 3.5.1);  mass = {massname} ===")
    print(f"  {'n':>5} {'R=ext/n':>12} {'sigma':>12} {'min cell':>12} "
          f"{'max cell':>12} {'empty':>8}")
    for n in args.cells:
        s = sigma_grid(dom["points"], dom["edges"], m_edge, n)
        flag = "" if s["R"] >= R0 else "  (R < R0!)"
        print(f"  {s['n']:>5} {s['R']:>12.4g} {s['sigma']:>12.4g} "
              f"{s['cmin']:>12.4g} {s['cmax']:>12.4g} "
              f"{s['n_empty']:>4}/{s['n_cells']}{flag}")
    print("  note: sigma at R < R0 is below the microstructure scale and only "
          "reflects\n        discretization; the assumption is stated for R >= R0.")

    # ---- mu (optional) ---------------------------------------------------
    if args.mu:
        print()
        print("=== mu  (connectivity, Lemma 3.6, Neumann; estimate) ===")
        print(f"  {'n':>5} {'R=ext/n':>12} {'mu':>12} {'cells eval':>12} "
              f"{'disc':>6} {'skip':>6}")
        for n in args.cells:
            m = estimate_mu(dom["points"], dom["edges"], m_edge, lengths, R0, n,
                            max_nodes=args.mu_max_nodes)
            print(f"  {m['n']:>5} {m['R']:>12.4g} {m['mu']:>12.4g} "
                  f"{m['n_eval']:>12} {m['n_disconnected']:>6} {m['n_skip']:>6}")
        print("  note: induced-subgraph estimate on the R0-enlarged cell (not the "
              "BFS-minimal\n        subgraph of the paper); disc = cells whose "
              "subgraph was disconnected.")


if __name__ == "__main__":
    main()
