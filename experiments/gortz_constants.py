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
      over the cells T of a regular n x n grid on the domain.  B_R(x) is the
      box of RADIUS R around x (side 2R), so n cells per axis correspond to
      R = 1/(2n): the paper's column R^-1 equals 2n.
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
cell B_R(x) from the generalized eigenproblem  Lbar u = lambda_2 Mbar u
(Neumann) on the connected component covering B_R(x) of the subgraph induced
on the R0-enlarged box, mu = max_x lambda_2^-1/2 / (2R).  See estimate_mu()
for the normalization and the relation to the paper's BFS-minimal subgraph.

Usage:
    python gortz_constants.py domain.geo.h5
    python gortz_constants.py domain.geo.h5 --cells 2 4 8 16 32
    python gortz_constants.py domain.geo.h5 --use-properties-mass
    python gortz_constants.py domain.geo.h5 --mu --cells 4 8 16
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
    """sigma(R) on a regular n x n tiling (box side 2R = extent/n) of the xy bbox.

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
    return dict(n=n, R=float(ext.max() / (2 * n)), sigma=sigma,
                cmin=float(pos.min()) if pos.size else 0.0,
                cmax=float(cell_mass.max()),
                mean=float(cell_mass.mean()), n_cells=n * n, n_empty=n_empty)


def estimate_mu(points, edges, m_edge, lengths, R0, n, max_nodes=200000):
    """Connectivity constant mu of Lemma 3.6, Neumann case, per cell.

    For each grid cell B_R(x) (a box of side 2R) we take the subgraph induced
    on the R0-enlarged box B_{R+R0}(x), split it into connected components,
    and use the component containing B_R(x)'s nodes.  Detached fragments that
    merely clip the enlarged box's boundary are irrelevant to the assumption
    and are discarded; if the CORE nodes themselves fall into more than one
    component, the connectivity assumption is genuinely violated at this cell
    (counted in n_violation, worst core coverage in min_coverage) and the
    dominant component is used.  The paper instead grows an edge-minimal
    connected subgraph by BFS; the full component can only be better
    connected, so this estimate is a slightly optimistic (lower) mu.

    On that component: weighted graph Laplacian Lbar (edge weight 1/|e|,
    operator L of Example 3.2 with gamma=1), mass Mbar lumped from the
    component's own edges, lambda_2 = second-smallest eigenvalue of
    Lbar u = lambda Mbar u.  Then  mu(x,R) = lambda_2^-1/2 / (2R)  and
    mu = max over cells.  Normalizing by the box side 2R reproduces the
    paper's Table 1; reading their inequality lambda_2^-1 <= mu^2 R^2 with R
    the box radius would double every value.  The average lambda_2^-1 is
    also returned for comparison with the paper's Figure 6.

    Cells whose component exceeds max_nodes are skipped (reported)."""
    from scipy.sparse.linalg import eigsh

    xy = points[:, :2]
    lo = xy.min(0)
    ext = xy.max(0) - lo
    ext = np.where(ext > 0, ext, 1.0)
    hx, hy = ext / n
    R = float(ext.max() / (2 * n))          # box radius: cells have side 2R

    w = 1.0 / lengths                       # Laplacian edge weights

    mu_max, worst_cell = 0.0, None
    mus, invlam2 = [], []
    n_skip, n_viol, n_eval = 0, 0, 0
    min_cov = 1.0
    for i in range(n):
        for j in range(n):
            x0, x1 = lo[0] + i * hx, lo[0] + (i + 1) * hx
            y0, y1 = lo[1] + j * hy, lo[1] + (j + 1) * hy
            # nodes of the R0-enlarged cell
            sel = ((xy[:, 0] >= x0 - R0) & (xy[:, 0] <= x1 + R0) &
                   (xy[:, 1] >= y0 - R0) & (xy[:, 1] <= y1 + R0))
            idx = np.nonzero(sel)[0]
            if idx.size < 3:
                continue

            remap = -np.ones(len(points), dtype=np.int64)
            remap[idx] = np.arange(idx.size)
            emask = sel[edges[:, 0]] & sel[edges[:, 1]]
            se = remap[edges[emask]]
            if se.size == 0:
                continue
            we = w[emask]
            me = m_edge[emask]

            A = sp.csr_matrix((we, (se[:, 0], se[:, 1])),
                              shape=(idx.size, idx.size))
            n_comp, labels = sp.csgraph.connected_components(A, directed=False)

            # component covering the core cell
            core = ((xy[idx, 0] >= x0) & (xy[idx, 0] <= x1) &
                    (xy[idx, 1] >= y0) & (xy[idx, 1] <= y1))
            if not core.any():
                continue
            core_labels = labels[core]
            dominant = np.bincount(core_labels).argmax()
            cov = float((core_labels == dominant).mean())
            if cov < 1.0:                   # core split across components:
                n_viol += 1                 # connectivity assumption violated
                min_cov = min(min_cov, cov)

            keep = labels == dominant
            n_sub = int(keep.sum())
            if n_sub > max_nodes:
                n_skip += 1
                continue
            if n_sub < 3:
                continue
            sub = -np.ones(idx.size, dtype=np.int64)
            sub[keep] = np.arange(n_sub)
            gmask = keep[se[:, 0]]          # edges stay within one component
            ge, gw, gm = sub[se[gmask]], we[gmask], me[gmask]

            # weighted Laplacian L = D - W, mass Mbar lumped from own edges
            W = sp.coo_matrix((gw, (ge[:, 0], ge[:, 1])), shape=(n_sub, n_sub))
            W = (W + W.T)
            L = sp.diags(np.asarray(W.sum(1)).ravel()) - W
            node_mass = np.zeros(n_sub)
            np.add.at(node_mass, ge[:, 0], gm / 2.0)
            np.add.at(node_mass, ge[:, 1], gm / 2.0)
            if node_mass.min() <= 0:
                continue                    # ill-defined Mbar
            M = sp.diags(node_mass)

            try:
                # smallest two eigenvalues of the SPD pencil; shift just below 0
                vals = eigsh(L.tocsc(), k=2, M=M.tocsc(),
                             sigma=-1e-12 * gw.sum(), which="LM",
                             return_eigenvectors=False)
                lam2 = float(np.sort(vals)[1])
            except Exception:
                continue
            if lam2 <= 0:
                continue

            n_eval += 1
            mu_cell = lam2 ** -0.5 / (2 * R)
            mus.append(mu_cell)
            invlam2.append(1.0 / lam2)
            if mu_cell > mu_max:
                mu_max, worst_cell = mu_cell, (i, j)

    return dict(n=n, R=R, mu=mu_max,
                mu_mean=float(np.mean(mus)) if mus else 0.0,
                avg_invlam2=float(np.mean(invlam2)) if invlam2 else 0.0,
                worst_cell=worst_cell, n_eval=n_eval,
                n_violation=n_viol, min_coverage=min_cov, n_skip=n_skip)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("domain", help="input .geo.h5 domain (make_geo2 format)")
    ap.add_argument("--cells", type=int, nargs="+", default=[2, 4, 8, 16, 32],
                    metavar="N",
                    help="cells per axis for the sigma sweep; the cells are the "
                         "boxes B_R(x) of radius R = extent/(2N), i.e. the "
                         "paper's R^-1 = 2N; default reproduces R^-1 = 4..64")
    ap.add_argument("--use-properties-mass", action="store_true",
                    help="use properties column 0 as edge mass instead of the "
                         "geometric edge length (default)")
    ap.add_argument("--mu", action="store_true",
                    help="also estimate the connectivity constant mu (slow)")
    ap.add_argument("--R0", type=float, default=None,
                    help="length scale R0 enlarging the cells in the mu "
                         "estimate; default: the computed R0 (paper uses 1/64)")
    ap.add_argument("--mu-max-nodes", type=int, default=200000,
                    help="skip cells whose connected component exceeds this "
                         "many nodes in the mu estimate (default 200000)")
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
    print(f"  {'n':>5} {'R=ext/2n':>12} {'R^-1':>8} {'sigma':>12} {'min cell':>12} "
          f"{'max cell':>12} {'empty':>8}")
    for n in args.cells:
        s = sigma_grid(dom["points"], dom["edges"], m_edge, n)
        flag = "" if s["R"] >= R0 else "  (R < R0!)"
        print(f"  {s['n']:>5} {s['R']:>12.4g} {1.0 / s['R']:>8.4g} {s['sigma']:>12.4g} "
              f"{s['cmin']:>12.4g} {s['cmax']:>12.4g} "
              f"{s['n_empty']:>4}/{s['n_cells']}{flag}")
    print("  note: sigma at R < R0 is below the microstructure scale and only "
          "reflects\n        discretization; the assumption is stated for R >= R0.")

    # ---- mu (optional) ---------------------------------------------------
    if args.mu:
        R0_mu = args.R0 if args.R0 is not None else R0
        print()
        print(f"=== mu  (connectivity, Lemma 3.6, Neumann; estimate)  "
              f"R0 = {R0_mu:g} ({'--R0' if args.R0 is not None else 'computed'}) ===")
        print(f"  {'n':>5} {'R=ext/2n':>12} {'R^-1':>8} {'mu':>10} {'mu_mean':>10} "
              f"{'avg lam2^-1':>12} {'eval':>6} {'viol':>6} {'skip':>6}")
        for n in args.cells:
            m = estimate_mu(dom["points"], dom["edges"], m_edge, lengths, R0_mu, n,
                            max_nodes=args.mu_max_nodes)
            print(f"  {m['n']:>5} {m['R']:>12.4g} {1.0 / m['R']:>8.4g} {m['mu']:>10.4g} "
                  f"{m['mu_mean']:>10.4g} {m['avg_invlam2']:>12.4g} "
                  f"{m['n_eval']:>6} {m['n_violation']:>6} {m['n_skip']:>6}")
        print("  note: lambda_2 on the connected component covering B_R(x) of the "
              "subgraph\n        induced on the R0-enlarged box; mu = lambda_2^-1/2 "
              "/ (2R), the side\n        normalization that reproduces the paper's "
              "Table 1.  viol = cells whose\n        core nodes span several "
              "components (connectivity violated there).")


if __name__ == "__main__":
    main()
