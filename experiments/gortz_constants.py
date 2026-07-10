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
    python gortz_constants.py domain.geo.h5 --mu --cells 4 8 16 --jobs 16
"""

import argparse
import os
import time
import sys

# force single-threaded BLAS/OpenMP before numpy loads: Arch's OpenBLAS is
# OpenMP-built and libgomp is not fork-safe -- once the parent has run any
# threaded kernel, eigsh in a fork()ed mu worker deadlocks.  Parallelism is
# over cells (--jobs), not inside the eigensolves, so nothing is lost.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

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
        mass = stiff = fiber_id = None
        if "properties" in g:
            P = g["properties"][:]
            mass = P[:, 0].astype(float)
            if P.shape[1] >= 7:  # [mass, EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2, ...]
                stiff = P[:, 1:7].astype(float)
            if P.shape[1] >= 16:
                fiber_id = P[:, 15].astype(np.int64)
        types_points = None
        if "types_points" in g:
            types_points = np.asarray(g["types_points"][:]).reshape(-1)
        size = np.array(g.attrs["size"]) if "size" in g.attrs else None
    tprint(f"loaded '{path}': {len(points)} nodes, {len(edges)} edges, "
           f"properties={'yes' if mass is not None else 'no'}")
    return dict(points=points, edges=edges, mass=mass, stiff=stiff,
                fiber_id=fiber_id, types_points=types_points, size=size)


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


def compute_lc(stiff, fiber_id=None):
    """Per-fiber Timoshenko bending length l_c = sqrt(EI/GA), the crossover
    span between shear- and bending-dominated transverse response.  The shear
    term GA |u' - r x t|^2 is the only coupling between the trace fields u and
    r; over a span R it penalizes a componentwise-interpolated (u, r) pair by
    ~ GA (amp/R)^2 * mass, while the low-energy bending mode it has to
    approximate only costs ~ EI (amp/R^2)^2 * mass, smaller by (l_c/R)^2.  A
    coarse space that interpolates u and r independently per component
    (net2as q1/pu) therefore loses the stable-decomposition property by the
    factor (R/l_c)^2: for R >> l_c the coarse level is inert and the method
    behaves one-level (kappa ~ R^-2, its ~ R^-1; measured on fiber2 in
    ne18-16..18).  Robust gortz-style behaviour needs R <~ l_c, so a regime
    window exists iff l_c >= R0.  stiff columns: [EA, kG_1A, kG_2A, G_xI_x,
    E_1I_1, E_2I_2]; l_c is a cross-section property, constant along a fiber,
    so percentiles are taken over unique fibers when fiber_id is present."""
    if fiber_id is not None:
        _, first = np.unique(fiber_id, return_index=True)
        stiff = stiff[first]
    ga1, ga2, ei1, ei2 = stiff[:, 1], stiff[:, 2], stiff[:, 4], stiff[:, 5]
    valid = (ga1 > 0) & (ga2 > 0) & (ei1 > 0) & (ei2 > 0)
    lc1 = np.sqrt(ei1[valid] / ga1[valid])
    lc2 = np.sqrt(ei2[valid] / ga2[valid])
    lc_min = np.minimum(lc1, lc2)
    pct = lambda v: {p: np.percentile(v, p) for p in (2, 50, 98)}
    return dict(lc1=pct(lc1), lc2=pct(lc2), lc_min=pct(lc_min),
                n_fibers=len(lc_min), n_invalid=int((~valid).sum()))


def compute_wave_speeds(mass, stiff, lengths):
    """Per-edge Timoshenko wave speeds c = sqrt(stiffness / (mass/length)):
    longitudinal c_l = sqrt(EA / mu) and shear c_s = sqrt(min(kG_1A, kG_2A) / mu)
    with mu = properties mass / edge length (the physical line density, NOT the
    paper's density-1 convention).  Material edges only (mass > 0): virtual weld
    edges are massless, respond quasi-statically and carry no waves of their
    own.  stiff columns: [EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2]."""
    mu = mass / lengths
    valid = mu > 0
    cl = np.sqrt(stiff[valid, 0] / mu[valid])
    cs = np.sqrt(np.minimum(stiff[valid, 1], stiff[valid, 2]) / mu[valid])
    pct = lambda v: {p: float(np.percentile(v, p)) for p in (2, 50, 98)}
    return dict(cl=pct(cl), cs=pct(cs),
                n_edges=int(valid.sum()), n_massless=int((~valid).sum()))


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


_MU = None   # read-only per-call state, inherited copy-on-write by fork()ed workers


def _mu_cell(ij):
    """Evaluate one grid cell (i, j) of the mu estimate (see estimate_mu).

    Reads the shared state in _MU.  Returns None for cells without a usable
    subgraph, else (i, j, coverage, skipped, mu_cell, 1/lambda_2) where the
    last two are None when the eigenproblem was not solved."""
    from scipy.sparse.linalg import eigsh

    g = _MU
    i, j = ij
    xy, R0 = g["xy"], g["R0"]
    x0, x1 = g["lo"][0] + i * g["hx"], g["lo"][0] + (i + 1) * g["hx"]
    y0, y1 = g["lo"][1] + j * g["hy"], g["lo"][1] + (j + 1) * g["hy"]

    # nodes of the R0-enlarged cell: x-band from the presorted x-coordinates
    # (two searchsorted calls instead of a full-array mask), y test on the band
    a = np.searchsorted(g["xs"], x0 - R0, side="left")
    b = np.searchsorted(g["xs"], x1 + R0, side="right")
    cand = g["xorder"][a:b]
    yc = xy[cand, 1]
    idx = cand[(yc >= y0 - R0) & (yc <= y1 + R0)]
    if idx.size < 3:
        return None
    idx.sort()

    sel = np.zeros(g["n_points"], dtype=bool)
    sel[idx] = True

    # candidate edges from the same x-band trick: an edge with both endpoints
    # selected has its min endpoint x inside [x0-R0, x1+R0]
    ea = np.searchsorted(g["exs"], x0 - R0, side="left")
    eb = np.searchsorted(g["exs"], x1 + R0, side="right")
    ec = g["edges_x"][ea:eb]
    emask = sel[ec[:, 0]] & sel[ec[:, 1]]
    if not emask.any():
        return None
    se = np.searchsorted(idx, ec[emask])    # global -> local node ids
    we = g["w_x"][ea:eb][emask]
    me = g["m_x"][ea:eb][emask]

    A = sp.csr_matrix((we, (se[:, 0], se[:, 1])), shape=(idx.size, idx.size))
    n_comp, labels = sp.csgraph.connected_components(A, directed=False)

    # component covering the core cell
    core = ((xy[idx, 0] >= x0) & (xy[idx, 0] <= x1) &
            (xy[idx, 1] >= y0) & (xy[idx, 1] <= y1))
    if not core.any():
        return None
    core_labels = labels[core]
    dominant = np.bincount(core_labels).argmax()
    cov = float((core_labels == dominant).mean())

    keep = labels == dominant
    n_sub = int(keep.sum())
    if n_sub > g["max_nodes"]:
        return (i, j, cov, True, None, None)
    if n_sub < 3:
        return (i, j, cov, False, None, None)
    sub = -np.ones(idx.size, dtype=np.int64)
    sub[keep] = np.arange(n_sub)
    gmask = keep[se[:, 0]]              # edges stay within one component
    ge, gw, gm = sub[se[gmask]], we[gmask], me[gmask]

    # weighted Laplacian L = D - W, mass Mbar lumped from own edges
    W = sp.coo_matrix((gw, (ge[:, 0], ge[:, 1])), shape=(n_sub, n_sub))
    W = (W + W.T)
    L = sp.diags(np.asarray(W.sum(1)).ravel()) - W
    node_mass = np.zeros(n_sub)
    np.add.at(node_mass, ge[:, 0], gm / 2.0)
    np.add.at(node_mass, ge[:, 1], gm / 2.0)
    if node_mass.min() <= 0:
        return (i, j, cov, False, None, None)   # ill-defined Mbar
    M = sp.diags(node_mass)

    try:
        # smallest two eigenvalues of the SPD pencil; shift just below 0
        vals = eigsh(L.tocsc(), k=2, M=M.tocsc(),
                     sigma=-1e-12 * gw.sum(), which="LM",
                     return_eigenvectors=False)
        lam2 = float(np.sort(vals)[1])
    except Exception:
        return (i, j, cov, False, None, None)
    if lam2 <= 0:
        return (i, j, cov, False, None, None)

    mu_cell = lam2 ** -0.5 / (2 * g["R"])
    return (i, j, cov, False, mu_cell, 1.0 / lam2)


def estimate_mu(points, edges, m_edge, lengths, R0, n, max_nodes=200000, jobs=1):
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

    Cells whose component exceeds max_nodes are skipped (reported).

    Cells are independent and evaluated by a fork()ed process pool of `jobs`
    workers (largest cells first, so a big straggler doesn't hold the tail);
    the shared arrays are inherited copy-on-write via the module global _MU."""
    global _MU

    xy = points[:, :2]
    lo = xy.min(0)
    ext = xy.max(0) - lo
    ext = np.where(ext > 0, ext, 1.0)
    hx, hy = ext / n
    R = float(ext.max() / (2 * n))          # box radius: cells have side 2R

    # presort nodes and edges by x so each cell extracts its x-band with two
    # searchsorted calls instead of scanning all points/edges
    xorder = np.argsort(xy[:, 0])
    exlo = np.minimum(xy[edges[:, 0], 0], xy[edges[:, 1], 0])
    eorder = np.argsort(exlo)
    _MU = dict(xy=xy, n_points=len(points), R0=R0, lo=lo, hx=hx, hy=hy, R=R,
               max_nodes=max_nodes, xorder=xorder, xs=xy[xorder, 0],
               exs=exlo[eorder], edges_x=edges[eorder],
               w_x=(1.0 / lengths)[eorder],     # Laplacian edge weights
               m_x=m_edge[eorder])

    # schedule the expensive cells first: node count of the R0-enlarged box
    # from a summed-area table of the per-cell node counts
    ix = np.clip(((xy[:, 0] - lo[0]) / ext[0] * n).astype(int), 0, n - 1)
    iy = np.clip(((xy[:, 1] - lo[1]) / ext[1] * n).astype(int), 0, n - 1)
    cnt = np.bincount(ix * n + iy, minlength=n * n).reshape(n, n)
    P = np.zeros((n + 1, n + 1))
    P[1:, 1:] = cnt.cumsum(0).cumsum(1)
    kx, ky = int(np.ceil(R0 / hx)), int(np.ceil(R0 / hy))

    def box_count(ij):
        i0, i1 = max(ij[0] - kx, 0), min(ij[0] + kx + 1, n)
        j0, j1 = max(ij[1] - ky, 0), min(ij[1] + ky + 1, n)
        return P[i1, j1] - P[i0, j1] - P[i1, j0] + P[i0, j0]

    cells = sorted(((i, j) for i in range(n) for j in range(n)),
                   key=box_count, reverse=True)

    if jobs > 1:
        import multiprocessing
        recs = []
        step = max(1, len(cells) // 10)
        with multiprocessing.get_context("fork").Pool(jobs) as pool:
            for k, rec in enumerate(pool.imap_unordered(_mu_cell, cells,
                                                        chunksize=1), 1):
                recs.append(rec)
    else:
        recs = [_mu_cell(c) for c in cells]

    # reduce in (i, j) order so the result is independent of completion order
    mu_max, worst_cell = 0.0, None
    mus, invlam2 = [], []
    n_skip, n_viol, n_eval = 0, 0, 0
    min_cov = 1.0
    for i, j, cov, skipped, mu_cell, il2 in sorted(r for r in recs
                                                   if r is not None):
        if cov < 1.0:                       # core split across components:
            n_viol += 1                     # connectivity assumption violated
            min_cov = min(min_cov, cov)
        if skipped:
            n_skip += 1
        if mu_cell is None:
            continue
        n_eval += 1
        mus.append(mu_cell)
        invlam2.append(il2)
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
    ap.add_argument("--jobs", type=int, default=os.cpu_count(),
                    help="worker processes for the mu estimate "
                         "(default: all cores)")
    ap.add_argument("--csv", metavar="PATH", default=None,
                    help="additionally write the per-n table (n, R, sigma, R/l_c, mu when "
                         "--mu) as CSV; scalars (extent, R0, l_c percentiles) go into "
                         "'# key = value' header comment lines (pandas: comment='#', "
                         "pgfplotstable: comment chars=#)")
    args = ap.parse_args()

    # per-n rows and scalar header for --csv, filled alongside the printed sections
    csv_rows = {n: dict(n=n) for n in args.cells}
    csv_scalars = {}

    dom = load_domain(args.domain)
    lengths = edge_lengths(dom["points"], dom["edges"])
    m_edge = edge_mass(dom, args.use_properties_mass, lengths)

    xy_ext = (dom["points"][:, :2].max(0) - dom["points"][:, :2].min(0))
    tprint(f"xy bounding box extent: {xy_ext[0]:.4g} x {xy_ext[1]:.4g}")
    csv_scalars["extent_x"] = xy_ext[0]
    csv_scalars["extent_y"] = xy_ext[1]

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
    csv_scalars["R0"] = R0

    # ---- l_c (Timoshenko bending length; not a constant of gortz.pdf) -----
    print()
    print("=== l_c  (Timoshenko bending length sqrt(EI/GA); componentwise "
          "coarse-space regime) ===")
    if dom["stiff"] is None:
        print("  no stiffness properties in file: solver defaults to unit "
              "stiffnesses => l_c = 1 (coordinate units)")
        lc = {p: 1.0 for p in (2, 50, 98)}
    else:
        lci = compute_lc(dom["stiff"], dom["fiber_id"])
        lc = lci["lc_min"]
        for nm, v in (("l_c1 (dir 1)", lci["lc1"]), ("l_c2 (dir 2)", lci["lc2"]),
                      ("min(l_c1,l_c2)", lci["lc_min"])):
            print(f"  {nm:<15}: p2={v[2]:.4g}  p50={v[50]:.4g}  p98={v[98]:.4g}"
                  f"   ({lci['n_fibers']} fibers"
                  + (f", {lci['n_invalid']} nonpos skipped" if lci["n_invalid"]
                     else "") + ")")
    ok = "yes" if lc[50] >= R0 else "NO"
    print(f"  regime window l_c >= R0: {ok}   (l_c p50 = {lc[50]:.4g}, "
          f"R0 = {R0:.4g})")
    for p in (2, 50, 98):
        csv_scalars[f"lc_min_p{p}"] = lc[p]
    print(f"  {'n':>5} {'R=ext/2n':>12} {'R/lc_p50':>10} {'(R/lc_p50)^2':>13} "
          f"{'R/lc_p2':>10} {'(R/lc_p2)^2':>12}")
    for n in args.cells:
        R = xy_ext.max() / (2 * n)
        print(f"  {n:>5} {R:>12.4g} {R / lc[50]:>10.4g} {(R / lc[50])**2:>13.4g} "
              f"{R / lc[2]:>10.4g} {(R / lc[2])**2:>12.4g}")
        csv_rows[n].update(R=R, Rinv=1.0 / R, R_over_lc_p50=R / lc[50],
                           R_over_lc_p50_sq=(R / lc[50])**2,
                           R_over_lc_p2=R / lc[2], R_over_lc_p2_sq=(R / lc[2])**2)
    print("  note: (R/l_c)^2 multiplies the stable-decomposition bound of a "
          "coarse space\n        that interpolates u and r componentwise "
          "(net2as q1/pu).  R >> l_c =>\n        coarse level inert, one-level "
          "behaviour (kappa ~ R^-2).  Flat\n        gortz-style curves need "
          "R <~ l_c on every subdomain.")

    # ---- wave speeds / time scale (not a constant of gortz.pdf) -----------
    print()
    print("=== wave speeds / domain time scale (dataset units, morgan: kg-um-s) ===")
    if dom["stiff"] is None or dom["mass"] is None:
        print("  no mass/stiffness properties in file: skipped")
    else:
        ws = compute_wave_speeds(dom["mass"], dom["stiff"], lengths)
        for nm, v in (("c_long = sqrt(EA/(m/l))", ws["cl"]),
                      ("c_shear = sqrt(kGA/(m/l))", ws["cs"])):
            print(f"  {nm:<26}: p2={v[2]:.4g}  p50={v[50]:.4g}  p98={v[98]:.4g}"
                  f"   ({ws['n_edges']} material edges"
                  + (f", {ws['n_massless']} massless skipped" if ws["n_massless"]
                     else "") + ")")
        cs, cl = ws["cs"][50], ws["cl"][50]
        t_cross = xy_ext.max() / cs
        # transverse fundamental of the clamped extent from the Timoshenko dispersion:
        # shear branch om_s = c_s k, bending branch om_b = sqrt(EI/mu) k^2 = c_s l_c k^2
        # (l_c = sqrt(EI/GA) from the section above -- so this IS sensitive to
        # --rescale-props on the rotation rigidities, unlike c_shear itself); the
        # compliances add in series, 1/om^2 = 1/om_s^2 + 1/om_b^2.  k = pi sqrt(2)/L
        # for the doubly clamped square.  Regime: shear iff l_c k >= 1.
        k1 = np.pi * np.sqrt(2.0) / xy_ext.max()
        om_s = cs * k1
        om_b = cs * lc[50] * k1 * k1
        om1 = (om_s**-2 + om_b**-2) ** -0.5
        T1 = 2 * np.pi / om1
        regime = "shear" if lc[50] * k1 >= 1 else "bending"
        print(f"  domain crossing extent/c_shear_p50   : {t_cross:.4g}")
        print(f"  fundamental (k = pi*sqrt2/extent)    : T_shear={2*np.pi/om_s:.4g}  "
              f"T_bend={2*np.pi/om_b:.4g}  -> T1 = {T1:.4g}  ({regime}-dominated, "
              f"l_c*k = {lc[50]*k1:.3g})")
        print("  note: ne18-22 measured T1 = 2.0e-6 on the reg-x1e6 quarter vs 1.7e-6 "
              "here.  In the\n        bending regime the domain-scale beam estimate "
              "ignores network truss action\n        (off-axis fibers carry transverse "
              "load axially), so T_bend is an UPPER\n        estimate; segment-scale "
              "bending in series with shear can be much faster.")
        csv_scalars["c_long_p50"] = cl
        csv_scalars["c_shear_p50"] = cs
        csv_scalars["t_cross"] = t_cross
        csv_scalars["T1_shear"] = 2 * np.pi / om_s
        csv_scalars["T1_bend"] = 2 * np.pi / om_b
        csv_scalars["T1_est"] = T1

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
        csv_rows[n].update(sigma=s["sigma"], cell_min=s["cmin"], cell_max=s["cmax"],
                           cells_empty=s["n_empty"], cells_total=s["n_cells"])
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
                            max_nodes=args.mu_max_nodes, jobs=args.jobs)
            print(f"  {m['n']:>5} {m['R']:>12.4g} {1.0 / m['R']:>8.4g} {m['mu']:>10.4g} "
                  f"{m['mu_mean']:>10.4g} {m['avg_invlam2']:>12.4g} "
                  f"{m['n_eval']:>6} {m['n_violation']:>6} {m['n_skip']:>6}")
            csv_rows[n].update(mu=m["mu"], mu_mean=m["mu_mean"],
                               avg_invlam2=m["avg_invlam2"], mu_eval=m["n_eval"],
                               mu_viol=m["n_violation"], mu_skip=m["n_skip"])
        print("  note: lambda_2 on the connected component covering B_R(x) of the "
              "subgraph\n        induced on the R0-enlarged box; mu = lambda_2^-1/2 "
              "/ (2R), the side\n        normalization that reproduces the paper's "
              "Table 1.  viol = cells whose\n        core nodes span several "
              "components (connectivity violated there).")

    # ---- csv (optional) ----------------------------------------------------
    if args.csv:
        import csv as csv_mod
        order = ["n", "R", "Rinv", "sigma", "cell_min", "cell_max", "cells_empty",
                 "cells_total", "R_over_lc_p50", "R_over_lc_p50_sq", "R_over_lc_p2",
                 "R_over_lc_p2_sq", "mu", "mu_mean", "avg_invlam2", "mu_eval",
                 "mu_viol", "mu_skip"]
        present = set().union(*(row.keys() for row in csv_rows.values()))
        fields = [c for c in order if c in present]
        with open(args.csv, "w", newline="") as fh:
            for k, v in csv_scalars.items():
                fh.write(f"# {k} = {v:.6g}\n")
            w = csv_mod.DictWriter(fh, fieldnames=fields, restval="")
            w.writeheader()
            for n in args.cells:
                w.writerow({k: v for k, v in csv_rows[n].items() if k in fields})
        tprint(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
