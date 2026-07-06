#!/usr/bin/env python

import numpy as np
import time
import argparse
import pandas
import h5py
import yaml
import scipy.sparse as sp
from scipy.spatial import cKDTree

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)

class Network:
  QUIRKS = ["morgan-2026-01-30"]

  def generate_grid(self, nx, ny):
    tprint(f"generating grid graph {nx} x {ny} on unit square")
    h_x = 1.0 / (nx - 1) if nx > 1 else 0.0
    h_y = 1.0 / (ny - 1) if ny > 1 else 0.0

    n_nodes = nx * ny
    nodes = np.zeros((n_nodes, 3))
    ii, jj = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    nodes[:, 0] = (jj * h_x).ravel()
    nodes[:, 1] = (ii * h_y).ravel()

    # horizontal edges: (i,j) -- (i,j+1), for j < nx-1
    h_src = (ii[:, :-1] * nx + jj[:, :-1]).ravel()
    h_dst = (ii[:, :-1] * nx + jj[:, :-1] + 1).ravel()

    # vertical edges: (i,j) -- (i+1,j), for i < ny-1
    v_src = (ii[:-1, :] * nx + jj[:-1, :]).ravel()
    v_dst = ((ii[:-1, :] + 1) * nx + jj[:-1, :]).ravel()

    edges = np.column_stack([
        np.concatenate([h_src, v_src]),
        np.concatenate([h_dst, v_dst]),
    ]).astype(np.int64)

    tprint("nodes", nodes.shape)
    tprint("edges", edges.shape)

    self.nodes = nodes
    self.edges = edges
    self.info = {"size": np.array([1.0, 1.0, 0.0])}
    self.edgeProps = None

  def generate_honeycomb(self, nx, ny, nz=1):
    """Honeycomb lattice with exactly nx*ny flat-top hexagons.

    Each row is nx hexagons in zigzag (centers alternate y by sqrt(3)/2);
    rows are stacked by sqrt(3) along y, sharing top/bottom edges.
    All boundary edges belong to a hexagon perimeter — no dangling stubs.
    """
    tprint(f"generating honeycomb graph {nx} x {ny} x {nz}")
    s3 = np.sqrt(3.0)

    # Hex center grid: row k, column m -> (1.5*m, (m%2)*s3/2 + k*s3)
    kk, mm = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    cx = 1.5 * mm.astype(float)
    cy = (mm % 2) * (s3 / 2) + kk.astype(float) * s3
    centers = np.stack([cx, cy], axis=-1)  # (ny, nx, 2)

    # 6 vertex offsets per flat-top hex, counterclockwise from the right vertex
    offs = np.array([
        [ 1.0,  0.0   ],
        [ 0.5,  s3/2  ],
        [-0.5,  s3/2  ],
        [-1.0,  0.0   ],
        [-0.5, -s3/2  ],
        [ 0.5, -s3/2  ],
    ])
    raw = (centers[:, :, None, :] + offs[None, None, :, :]).reshape(-1, 2)

    # Dedupe shared vertices by rounding to 1e-6 precision
    key = np.round(raw * 1e6).astype(np.int64)
    _, inv, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    nodes2 = np.zeros((counts.shape[0], 2))
    np.add.at(nodes2, inv, raw)
    nodes2 /= counts[:, None]
    vid = inv.reshape(ny, nx, 6)

    # 6 perimeter edges per hex, dedupe shared edges between adjacent hexes
    edges_raw = np.vstack([
      np.stack([vid[:, :, v].ravel(), vid[:, :, (v + 1) % 6].ravel()], axis=1)
      for v in range(6)
    ])
    edges2d = np.unique(np.sort(edges_raw, axis=1), axis=0).astype(np.int64)

    # Rescale 2D so x extent = 1
    nodes2 -= nodes2.min(axis=0)
    scale = 1.0 / nodes2[:, 0].max()
    nodes2 *= scale

    # Layer spacing: typical bond length in-plane after scaling.
    # All in-plane edges have length 1*scale (= s in original units), so dz = scale.
    dz = scale

    n_per_layer = nodes2.shape[0]
    n_nodes = n_per_layer * nz
    nodes = np.zeros((n_nodes, 3))
    layer_ids = np.arange(nz)
    # tile xy across layers
    nodes[:, :2] = np.tile(nodes2, (nz, 1))
    nodes[:, 2]  = np.repeat(layer_ids * dz, n_per_layer)

    # In-plane edges, replicated per layer with offset
    offsets = (np.arange(nz) * n_per_layer)[:, None, None]   # (nz,1,1)
    in_plane = edges2d[None, :, :] + offsets                 # (nz, n_e2d, 2)
    in_plane = in_plane.reshape(-1, 2)

    # Vertical edges between consecutive layers
    base = np.arange(n_per_layer)
    if nz > 1:
      base = np.arange(n_per_layer)
      v_src = np.concatenate([base + k*n_per_layer     for k in range(nz-1)])
      v_dst = np.concatenate([base + (k+1)*n_per_layer for k in range(nz-1)])
      vert = np.column_stack([v_src, v_dst])
    else:
      vert = np.empty((0, 2), dtype=np.int64)
    edges = np.vstack([in_plane, vert]).astype(np.int64)

    tprint("nodes", nodes.shape)
    tprint("edges", edges.shape)
    self.nodes = nodes
    self.edges = edges
    size = nodes.max(axis=0) - nodes.min(axis=0)
    self.info = {"size": size}
    self.edgeProps = None

  def generate_mikado(self, mass, r=0.05, seed=0, min_edge=None):
    """Random mikado-style fiber network on the unit square (gortz.pdf, sec. 6.1).

    Fibers of fixed length r are placed with midpoints uniform in
    [-r/2, 1+r/2]^2 and uniformly random rotation, clipped to the unit
    square, until the total fiber length reaches `mass`. Every pairwise
    fiber intersection becomes a node splitting both fibers, nodes closer
    than min_edge (default r*1e-4) are merged, and only the largest
    connected component is kept, so the result is a single connected graph.
    """
    tprint(f"generating mikado graph: mass={mass:g}, fiber length r={r:g}, seed={seed}")

    # continuum percolation of 2D sticks: two isotropic sticks of length r cross
    # iff their midpoint offset lies in a parallelogram of area r^2*sin(theta),
    # so the mean crossings per stick is k = n * <r^2 sin> = (2/pi)*n*r^2 with
    # stick density n = mass/r per unit area. A giant (domain-spanning) component
    # emerges above the numerically known threshold n*r^2 = mass*r ~ 5.64
    # (i.e. k ~ 3.59), Mertens & Moore, Phys. Rev. E 86, 061109 (2012).
    density = mass * r
    k_mean = 2.0 / np.pi * density
    tprint(f"stick density n*r^2 = mass*r = {density:.3g} = {density / 5.6373:.2g} x threshold 5.64, "
           f"~{k_mean:.3g} crossings per fiber")
    if density < 5.6373:
      tprint("  below percolation threshold: NO giant component expected, "
             "the largest component will only be a small local cluster")
    else:
      tprint("  above percolation threshold: giant component expected")

    rng = np.random.default_rng(seed)

    # place fibers (batched) until the clipped total length reaches `mass`
    seg_a, seg_b, seg_len = [], [], []
    total = 0.0
    while total < mass:
      n_batch = max(1024, int(1.2 * (mass - total) / r))
      mid = rng.uniform(-0.5 * r, 1.0 + 0.5 * r, size=(n_batch, 2))
      ang = rng.uniform(0.0, np.pi, size=n_batch)
      half = 0.5 * r * np.column_stack([np.cos(ang), np.sin(ang)])
      p0, p1 = mid - half, mid + half

      # Liang-Barsky clip to [0,1]^2: keep param t in [t0,t1] with p*t <= q per side
      t0, t1 = np.zeros(n_batch), np.ones(n_batch)
      keep = np.ones(n_batch, dtype=bool)
      for dim in range(2):
        d = p1[:, dim] - p0[:, dim]
        for p, q in ((-d, p0[:, dim]), (d, 1.0 - p0[:, dim])):
          with np.errstate(divide='ignore', invalid='ignore'):
            tq = q / p
          t0 = np.where(p < 0, np.maximum(t0, tq), t0)
          t1 = np.where(p > 0, np.minimum(t1, tq), t1)
          keep &= ~((p == 0) & (q < 0))
      keep &= t0 < t1
      a = (p0 + t0[:, None] * (p1 - p0))[keep]
      b = (p0 + t1[:, None] * (p1 - p0))[keep]
      l = np.linalg.norm(b - a, axis=1)
      seg_a.append(a); seg_b.append(b); seg_len.append(l)
      total += l.sum()

    a, b = np.vstack(seg_a), np.vstack(seg_b)
    lengths = np.concatenate(seg_len)
    n_fib = min(np.searchsorted(np.cumsum(lengths), mass) + 1, len(lengths))
    a, b = a[:n_fib], b[:n_fib]
    tprint(f"placed {n_fib} fibers, total length {lengths[:n_fib].sum():.6g}")

    # candidate pairs: fibers of length <= r can only meet if midpoints are within r
    tree = cKDTree(0.5 * (a + b))
    cand = tree.query_pairs(r, output_type='ndarray')
    i, j = cand[:, 0], cand[:, 1]
    di, dj, w = b[i] - a[i], b[j] - a[j], a[j] - a[i]
    cross2d = lambda u, v: u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    denom = cross2d(di, dj)
    ok = np.abs(denom) > 1e-12  # exactly parallel fibers (probability zero) are skipped
    with np.errstate(divide='ignore', invalid='ignore'):
      s = cross2d(w, dj) / denom
      t = cross2d(w, di) / denom
    ok &= (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
    i, j, s, t = i[ok], j[ok], s[ok], t[ok]
    pts = a[i] + s[:, None] * (b[i] - a[i])
    n_x = len(pts)
    tprint(f"{n_x} fiber-fiber intersections ({len(cand)} candidate pairs)")

    # nodes: fiber endpoints then intersection points; split each fiber into
    # edges between consecutive parameters along it
    raw_nodes = np.vstack([a, b, pts])
    fiber = np.concatenate([np.arange(n_fib), np.arange(n_fib), i, j])
    param = np.concatenate([np.zeros(n_fib), np.ones(n_fib), s, t])
    node  = np.concatenate([np.arange(2 * n_fib),
                            2 * n_fib + np.arange(n_x),
                            2 * n_fib + np.arange(n_x)])
    order = np.lexsort((param, fiber))
    fo, no = fiber[order], node[order]
    adj = fo[:-1] == fo[1:]
    edges = np.column_stack([no[:-1][adj], no[1:][adj]]).astype(np.int64)

    # merge close nodes, setting a lower bound on edge lengths; iterate since
    # merged centroids can again end up closer than the tolerance
    merge_tol = r * 1e-4 if min_edge is None else min_edge
    n_raw = raw_nodes.shape[0]
    merged = raw_nodes
    while True:
      close = cKDTree(merged).query_pairs(merge_tol, output_type='ndarray')
      if len(close) == 0:
        break
      n_cur = merged.shape[0]
      g = sp.csr_matrix((np.ones(len(close)), (close[:, 0], close[:, 1])), shape=(n_cur, n_cur))
      n_groups, labels = sp.csgraph.connected_components(g, directed=False)
      centroids = np.zeros((n_groups, 2))
      np.add.at(centroids, labels, merged)
      centroids /= np.bincount(labels)[:, None]
      edges = labels[edges]
      merged = centroids
    n_merged = merged.shape[0]
    tprint(f"merged {n_raw - n_merged} nodes closer than {merge_tol:.1e}")

    # drop self-loops and duplicate edges
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges[edges[:, 0] != edges[:, 1]], axis=0)

    # keep only the largest connected component
    A = sp.csr_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n_merged, n_merged))
    n_comp, comp = sp.csgraph.connected_components(A, directed=False)
    keep_node = comp == np.bincount(comp).argmax()
    remap = np.full(n_merged, -1, dtype=np.int64)
    remap[keep_node] = np.arange(keep_node.sum())
    edges = remap[edges[keep_node[edges[:, 0]]]]
    tprint(f"{n_comp} components, keeping largest with {keep_node.sum()} nodes")

    nodes = np.zeros((keep_node.sum(), 3))
    nodes[:, :2] = merged[keep_node]

    final_mass = np.linalg.norm(nodes[edges[:, 0]] - nodes[edges[:, 1]], axis=1).sum()
    tprint(f"final mass (total edge length): {final_mass:.6g}")
    tprint("nodes", nodes.shape)
    tprint("edges", edges.shape)
    self.nodes = nodes
    self.edges = edges
    self.info = {"size": np.array([1.0, 1.0, 0.0])}
    self.edgeProps = None

  def generate_synthetic_properties(self, width=None):
    """Build the 17-column edgeProps array for a synthetic network.

    density = 1 (mass = length), all stiffnesses = 1,
    n_1 = (0,0,1) (or fallback if tangent is vertical),
    n_2 = tangent x n_1 normalized, widths constant,
    fiber_id = 0..n_edges-1, fiber_edge_id = 0.

    If width is None, picks 0.1 * mean(edge length).
    """
    nodes = self.nodes
    edges = self.edges
    n_edges = edges.shape[0]
    tprint(f"generating synthetic properties for {n_edges} edges")

    p1 = nodes[edges[:, 0]]
    p2 = nodes[edges[:, 1]]
    tangent = p2 - p1
    lengths = np.linalg.norm(tangent, axis=1)
    t_hat = tangent / lengths[:, None]

    # n_1 = (0,0,1), projected orthogonal to tangent.
    # If the tangent is (nearly) parallel to z, fall back to (1,0,0).
    z = np.array([0.0, 0.0, 1.0])
    cos_tz = t_hat @ z
    near_z = np.abs(cos_tz) > 0.99
    n1 = np.tile(z, (n_edges, 1))
    n1[near_z] = np.array([1.0, 0.0, 0.0])
    # remove tangent component, renormalize
    n1 -= np.einsum('ij,ij->i', n1, t_hat)[:, None] * t_hat
    n1 /= np.linalg.norm(n1, axis=1, keepdims=True)
    # n_2 = tangent x n_1
    n2 = np.cross(t_hat, n1)
    n2 /= np.linalg.norm(n2, axis=1, keepdims=True)

    mass = lengths

    EA   = np.ones(n_edges)
    kG1A = np.ones(n_edges)
    kG2A = np.ones(n_edges)
    GxIx = np.ones(n_edges)
    E1I1 = np.ones(n_edges)
    E2I2 = np.ones(n_edges)

    if width is None:
      width = 0.1 * lengths.mean()
    tprint(f"  using width = {width:.3e}")
    width1 = np.full(n_edges, width)
    width2 = np.full(n_edges, width)

    fiber_id      = np.arange(n_edges, dtype=np.float64)
    fiber_edge_id = np.zeros(n_edges)

    self.edgeProps = np.column_stack([
      mass,
      EA, kG1A, kG2A,
      GxIx, E1I1, E2I2,
      n1, n2,
      width1, width2,
      fiber_id, fiber_edge_id,
    ])
    tprint("edgeProps", self.edgeProps.shape)

  def read_morgan(self, path, rescale_props=None, quirk=None):
    tprint("reading nodes")
    nodes   = pandas.read_csv(path + "/nodes.csv")
    nodes   = nodes.to_numpy()[:,1:]

    tprint("reading edges")
    edges   = pandas.read_csv(path + '/edges.csv')
    edges   = edges.to_numpy()[:,1:]

    tprint("reading edgeProps")
    edgeProps   = pandas.read_csv(path + '/edgeProperties.csv')
    edgeProps   = edgeProps.to_numpy()[:,1:]

    if rescale_props is not None:
      n_props = edgeProps.shape[-1]
      print(rescale_props)
      edgeProps *= rescale_props

    if quirk == "morgan-2026-01-30":
      n_edges = edges.shape[0]
      fiber_ids = np.arange(n_edges)
      fiber_edge_ids = np.zeros(n_edges)
      widths = np.ones(n_edges) * 10 # 10um default width
      edgeProps = np.column_stack([
        edgeProps, widths, widths, fiber_ids, fiber_edge_ids,
      ])

    tprint("nodes", nodes.shape)
    tprint("edges", edges.shape)
    tprint("edgeProps", edgeProps.shape)

    try:
      with open(path + '/info.txt') as f:
        info = yaml.safe_load(f)
    except FileNotFoundError:
      info = {}

    try:
      with open(path + '/units.txt') as f:
        info["units"] = f.read()
    except FileNotFoundError:
      pass

    mins = nodes.min(axis=0)
    maxs = nodes.max(axis=0)
    dims = maxs - mins
    info["size"] = dims

    self.nodes = nodes
    self.edges = edges
    self.edgeProps = edgeProps
    self.info = info

    # clamp to >= 1e-10
    #self.edgeProps = np.clip(self.edgeProps, 1e-10, 1e+10)

  def verify_nonzero(self):
    """Verify that material properties are nonzero where required.
    Reports per-column zero/near-zero counts and degenerate normal vectors.
    """
    if self.edgeProps is None:
      tprint("verify_nonzero: no edgeProps loaded, skipping")
      return

    eps = 1e-30
    props = self.edgeProps
    n = props.shape[0]
    tprint(f"verify_nonzero: checking {n} fibers")

    # column groups: (indices, label)
    groups = [
      ([0],              "mass"),
      ([1, 2, 3],        "displacement stiffness (EA, kG_1A, kG_2A)"),
      ([4, 5, 6],        "rotation stiffness (G_xI_x, E_1I_1, E_2I_2)"),
      ([13, 14],         "widths (width1, width2)"),
    ]

    all_ok = True
    for cols, label in groups:
      for c in cols:
        col = props[:, c]
        n_zero  = (col == 0).sum()
        n_small = ((np.abs(col) < eps) & (col != 0)).sum()
        n_neg   = (col < 0).sum()
        if n_zero or n_small or n_neg:
          all_ok = False
          tprint(f"  col {c:2d} ({label}): "
                 f"zero={n_zero} subnormal={n_small} negative={n_neg} "
                 f"min={col.min():.3e} max={col.max():.3e}")
        else:
          tprint(f"  col {c:2d} ({label}): ok "
                 f"min={col.min():.3e} max={col.max():.3e}")

    # normal vector lengths
    n1 = props[:, 7:10]
    n2 = props[:, 10:13]
    len1 = np.linalg.norm(n1, axis=1)
    len2 = np.linalg.norm(n2, axis=1)
    for vec_name, lens in [("normal 1", len1), ("normal 2", len2)]:
      n_zero  = (lens < eps).sum()
      n_nonunit = (np.abs(lens - 1.0) > 1e-6).sum()
      if n_zero or n_nonunit:
        all_ok = False
        tprint(f"  {vec_name} length: zero={n_zero} non-unit={n_nonunit} "
               f"min={lens.min():.3e} max={lens.max():.3e}")
      else:
        tprint(f"  {vec_name} length: ok (all unit)")

    fiber_id = props[:, 15].astype(np.int64)
    virtual = (fiber_id == -1)
    real    = ~virtual

    n_virt_zero_mass = (virtual & (props[:, 0] == 0)).sum()
    n_real_zero_mass = (real    & (props[:, 0] == 0)).sum()
    n_virt_total     = virtual.sum()
    n_real_total     = real.sum()

    tprint(f"virtual fibers (id=-1): {n_virt_total}, of which zero-mass: {n_virt_zero_mass}")
    tprint(f"real fibers:            {n_real_total}, of which zero-mass: {n_real_zero_mass}")

    # orthogonality of n1 and n2 (cheap bonus check)
    dots = np.einsum('ij,ij->i', n1, n2)
    n_nonorth = (np.abs(dots) > 1e-6).sum()
    if n_nonorth:
      all_ok = False
      tprint(f"  normal1 · normal2: non-orthogonal pairs={n_nonorth} "
               f"max|dot|={np.abs(dots).max():.3e}")
    else:
      tprint(f"  normal1 · normal2: ok (all orthogonal)")

    if all_ok:
      tprint("verify_nonzero: all checks passed")
    else:
      tprint("verify_nonzero: FAILED — see above")


  def prop_cutoff(self, percent):
    """Floor each stiffness component (edgeProps cols 1..6) at its `percent`-th percentile,
    shrinking the coefficient contrast (max/min) while raising at most ~`percent`% of fibers per
    component. For each component prints a per-decade log histogram, a few reference percentile
    thresholds for context, and the applied floor. Mutates edgeProps in place.
    """
    if self.edgeProps is None:
      tprint("prop_cutoff: no edgeProps loaded, skipping")
      return

    labels = ["mass", "EA", "kG_1A", "kG_2A", "G_xI_x", "E_1I_1", "E_2I_2"]
    props = self.edgeProps
    n = props.shape[0]
    if n == 0:
      tprint("prop_cutoff: 0 fibers, skipping")
      return
    ref = [0.5, 1, 2, 5]
    tprint(f"prop_cutoff: flooring cols 1..6 at the {percent:g}-th percentile, {n} fibers")

    for c in range(1, 7):
      col = props[:, c]
      cmin, cmax = col.min(), col.max()
      contrast = cmax / cmin if cmin > 0 else np.inf
      tprint(f"  col {c} ({labels[c]}): min={cmin:.3e} max={cmax:.3e} "
             f"contrast={contrast:.3e} nonpos={(col <= 0).sum()}")

      pos = col[col > 0]
      if pos.size:
        lo_e = int(np.floor(np.log10(pos.min())))
        hi_e = int(np.ceil(np.log10(pos.max())))
        edges = 10.0 ** np.arange(lo_e, hi_e + 1)
        hist, _ = np.histogram(pos, bins=edges)
        bars = " ".join(f"1e{lo_e+i:+03d}:{hist[i]}" for i in range(len(hist)))
        tprint(f"    log-hist (#fibers per decade): {bars}")

      refs = " ".join(f"p{p:g}={np.percentile(col, p):.3e}" for p in ref)
      tprint(f"    ref thresholds: {refs}")

      # apply the floor at the requested percentile
      t = np.percentile(col, percent)
      n_below = int((col < t).sum())
      props[:, c] = np.clip(col, t, None)
      contrast_after = cmax / t if t > 0 else np.inf
      tprint(f"    floored at p{percent:g}={t:.3e}: raised {n_below} fibers "
             f"({100*n_below/n:.3g}%), contrast {contrast:.3e} -> {contrast_after:.3e}")


  def node_edge_dedupe(self, merge_tol):
    nodes = self.nodes
    edges = self.edges
    edgeProps = self.edgeProps

    n_nodes = nodes.shape[0]
    n_edges = edges.shape[0]
    n_edgeProps = edgeProps.shape[0]
    edgeProps_dim = edgeProps.shape[1]

    _, inv, counts = np.unique(nodes, axis=0, return_inverse=True, return_counts=True)
    n_dup = (counts > 1).sum()
    tprint(f"nodes: bit exact duplicates: {n_dup}")

    canon = np.sort(edges, axis=1)
    unique_edges, edge_inv, edge_counts = np.unique(canon, axis=0, return_inverse=True, return_counts=True)
    n_dup_edges = (edge_counts > 1).sum()
    tprint(f"edges: duplicates: {n_dup_edges}, self-loops: {(edges[:,0] == edges[:,1]).sum()}")

    tprint("build KD tree")
    tree = cKDTree(nodes)
    d, _ = tree.query(nodes, k=2)   # k=1 is self → distance 0
    nn = d[:, 1]                     # nearest non-self distance

    tprint("nearest neighbor distance")
    tprint("percentile: distance")
    for p in [0.01, 0.1, 1, 5, 50, 95, 99]:
        tprint(f"  {p:>4}: {np.percentile(nn, p):.3e}")
    tprint(f"min={nn.min():.3e} max={nn.max():.3e}")
    tprint("<= distance: count")
    for tol in [1e-12, 1e-9, 1e-6, 1e-3, 1e-2, 1e-1, 1e0]:
        tprint(f"  {tol:.0e}: {(nn < tol).sum()}")

    close_pairs = tree.query_pairs(merge_tol, output_type='ndarray') # (n_pairs, 2), i < j
    tprint(f"merging {len(close_pairs)} pairs closer than {merge_tol:.3e}")

    g = sp.csr_matrix((np.ones(len(close_pairs)), (close_pairs[:,0], close_pairs[:,1])), shape=(n_nodes, n_nodes))
    n_groups, labels = sp.csgraph.connected_components(g, directed=False)

    new_nodes = np.zeros((n_groups, 3))
    np.add.at(new_nodes, labels, nodes)
    new_nodes /= np.bincount(labels)[:, None]

    edges = labels[edges]
    nodes = new_nodes
    n_nodes = n_groups

    tprint(f"after merge")
    tprint("nodes", nodes.shape)
    tprint("edges", edges.shape)
    tprint("edgeProps", edgeProps.shape)

    canon = np.sort(edges, axis=1)
    unique_edges, idx, edge_inv, edge_counts = np.unique(canon, axis=0, return_index=True, return_inverse=True, return_counts=True)
    n_dup_edges = (edge_counts > 1).sum()
    n_loops = (unique_edges[:,0] == unique_edges[:,1]).sum()
    tprint(f"edges: duplicates: {n_dup_edges}, self-loops: {n_loops}")

    # keep one representative per unique edge, drop self-loops
    keep = unique_edges[:,0] != unique_edges[:,1]
    edges     = unique_edges[keep]
    edgeProps = edgeProps[idx][keep]
    n_edges   = len(edges)
    n_edgeProps = edgeProps.shape[0]
    tprint(f"after dedup: {n_edges} edges")

    assert(n_edges == n_edgeProps)

    self.nodes = nodes
    self.edges = edges
    self.edgeProps = edgeProps


  def compute_types(self, tol):
    nodes = self.nodes
    edges = self.edges

    mins = nodes.min(axis=0)
    maxs = nodes.max(axis=0)
    dims = maxs - mins

    # small circle of nodes around the (xy) domain center
    center = (mins[:2] + maxs[:2]) / 2
    radius = args.center_radius * max(dims[0], dims[1])
    dist = np.linalg.norm(nodes[:, :2] - center, axis=1)

    sides = {
      "xmin": nodes[:,0] - mins[0] <= tol * dims[0],
      "xmax": maxs[0] - nodes[:,0] <= tol * dims[0],
      "ymin": nodes[:,1] - mins[1] <= tol * dims[1],
      "ymax": maxs[1] - nodes[:,1] <= tol * dims[1],
      "center": dist <= radius,
    }
    tprint(f"center circle: c={center}, r={radius:.3e}, {sides['center'].sum()} nodes")

    dir_side = np.array([sides[x.split('=')[0]]  for x in args.dirichlet])
    dir_desc = np.array([int(x.split('=')[1], 0) for x in args.dirichlet], dtype=np.int32)
    self.types_points = np.bitwise_or.reduce(dir_side * dir_desc[:, None], axis=0).astype(np.int32)

    count_dir = (self.types_points != 0).sum()
    frac_dir = count_dir / len(self.types_points)
    tprint("count dir", count_dir)
    tprint(f"frac dir {frac_dir:.5e}")

    self.types_faces = self.types_points[edges].astype(np.int32)


  def drop_floating_and_small_components(self, min_comp_size):
    nodes, edges, edgeProps, types_points, types_faces = self.nodes, self.edges, self.edgeProps, self.types_points, self.types_faces

    n_nodes = nodes.shape[0]
    n_edgeProps = edgeProps.shape[0]

    A = sp.csr_matrix((np.ones(len(edges)), (edges[:,0], edges[:,1])), shape=(n_nodes, n_nodes))
    n_comp, labels = sp.csgraph.connected_components(A, directed=False)
    free = sum(1 for c in range(n_comp) if types_points[labels == c].sum() == 0)
    tprint(f"{n_comp} components, {free} without any Dirichlet node")

    sizes = np.bincount(labels)
    order = np.argsort(sizes)[::-1]
    tprint(f"component sizes: {sizes[order].tolist()}")

    free_mask = np.array([types_points[labels == c].sum() == 0 for c in range(n_comp)])
    free_sizes = sizes[free_mask]
    if len(free_sizes):
        tprint(f"floating: count={len(free_sizes)} total_nodes={free_sizes.sum()} "
               f"min={free_sizes.min()} max={free_sizes.max()} mean={free_sizes.mean():.1f}")

    keep_comp = sizes >= min_comp_size  # boolean per component
    # drop floating components (no Dirichlet node) regardless of size:
    keep_comp &= ~free_mask
    keep_node = keep_comp[labels]
    tprint(f"keeping {keep_comp.sum()}/{n_comp} components, {keep_node.sum()}/{n_nodes} nodes")

    # remap: old node index -> new node index, -1 for dropped
    remap = np.full(n_nodes, -1, dtype=edges.dtype)
    remap[keep_node] = np.arange(keep_node.sum())

    self.nodes        = nodes[keep_node]
    self.types_points = types_points[keep_node]
    n_nodes      = len(nodes)

    keep_edge = keep_node[edges[:,0]] & keep_node[edges[:,1]]
    self.edges     = remap[edges[keep_edge]]
    self.edgeProps = edgeProps[keep_edge]
    self.types_faces  = types_faces[keep_edge]
    n_edges   = len(edges)
    tprint(f"after pruning: {n_nodes} nodes, {n_edges} edges")


  def write_h5(self, out, no_props=False):
    tprint(f"writing h5 file to '{out}'")
    with h5py.File(out, "w") as f:
      g = f.create_group("domain")
      g.create_dataset("points", data=self.nodes, compression="gzip")
      g.create_dataset("edges", data=self.edges, compression="gzip")
      if hasattr(self, "edgeProps") and self.edgeProps is not None and not no_props:
        g.create_dataset("properties", data=self.edgeProps, compression="gzip")
      g.create_dataset("types_points", data=self.types_points, compression="gzip")
      g.create_dataset("types_faces", data=self.types_faces, compression="gzip")
      for k, v in self.info.items():
        g.attrs[k] = v


  def clamp_xy(self, x, y):
    nodes, edges, edgeProps = self.nodes, self.edges, self.edgeProps
    n_nodes_old, n_edges_old = nodes.shape[0], edges.shape[0]

    xmin = 0
    ymin = 0
    xmax, ymax = self.info["size"][:2] * np.array([x, y])

    inside = ((nodes[:,0] >= xmin) & (nodes[:,0] <= xmax) &
              (nodes[:,1] >= ymin) & (nodes[:,1] <= ymax))

    both_in    = inside[edges[:,0]] & inside[edges[:,1]]
    both_out   = (~inside[edges[:,0]]) & (~inside[edges[:,1]])
    partial    = ~(both_in | both_out)

    tprint(f"clamp xy to [{xmin},{xmax}] x [{ymin},{ymax}]")
    tprint(f"  edges fully outside:    {both_out.sum()}")
    tprint(f"  edges partially outside: {partial.sum()}")
    tprint(f"  edges kept:             {both_in.sum()} / {n_edges_old}")

    edges     = edges[both_in]
    edgeProps = edgeProps[both_in]

    # drop now-unused nodes
    used = np.zeros(n_nodes_old, dtype=bool)
    used[edges.ravel()] = True
    remap = np.full(n_nodes_old, -1, dtype=edges.dtype)
    remap[used] = np.arange(used.sum())

    self.nodes     = nodes[used]
    self.edges     = remap[edges]
    self.edgeProps = edgeProps
    tprint(f"  nodes kept: {used.sum()} / {n_nodes_old}")
    mins = self.nodes.min(axis=0)
    maxs = self.nodes.max(axis=0)
    dims = maxs - mins
    self.info["size"] = dims


  def write_vtkhdf_view(self, out):
    """Add VTKHDF view to the .geo.h5 file: virtual Connectivity over /domain/edges,
    plus real Offsets, Types, and NumberOf* datasets.
    """
    path = out
    tprint(f"adding VTKHDF view to '{path}'")
    with h5py.File(path, "a") as f:
      if "VTKHDF" in f:
          del f["VTKHDF"]

      n_points = f["domain/points"].shape[0]
      n_cells  = f["domain/edges"].shape[0]
      n_conn   = 2 * n_cells

      root = f.create_group("VTKHDF")
      root.attrs.create("Version", [2, 0], dtype="int64")
      root.attrs.create("Type", np.bytes_("UnstructuredGrid"))

      root["Points"] = h5py.SoftLink("/domain/points")

      layout = h5py.VirtualLayout(shape=(n_conn,), dtype=f["domain/edges"].dtype)
      layout[...] = h5py.VirtualSource(path, "domain/edges",
                                       shape=f["domain/edges"].shape)
      root.create_virtual_dataset("Connectivity", layout)

      root.create_dataset(
          "Offsets",
          data=np.arange(0, n_conn + 2, 2, dtype=np.int64),
          compression="gzip",
      )

      root.create_dataset(
        "Types",
        data=np.full(n_cells, 3, dtype=np.uint8),
        compression="gzip",
      )

      root.create_dataset("NumberOfPoints",          data=np.array([n_points], dtype=np.int64))
      root.create_dataset("NumberOfCells",           data=np.array([n_cells],  dtype=np.int64))
      root.create_dataset("NumberOfConnectivityIds", data=np.array([n_conn],   dtype=np.int64))

      if "domain/types_points" in f:
        pd = root.create_group("PointData")
        pd["types_points"] = h5py.SoftLink("/domain/types_points")
      if "domain/properties" in f:
        cd = root.create_group("CellData")
        cd["properties"] = h5py.SoftLink("/domain/properties")


  def rescale_bbox(self):
    mins = self.nodes.min(axis=0)
    maxs = self.nodes.max(axis=0)
    dims = maxs - mins
    scale = 1.0 / max(dims[0], dims[1])
    tprint(f"rescaling by {scale:.3e} (bbox was {dims})")
    self.nodes = (self.nodes - mins) * scale
    self.info["size"] = (maxs - mins) * scale

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
  parser.add_argument("-i", "--input", help="input", default=".")
  parser.add_argument("-o", "--output", help="output", default="graph")
  parser.add_argument("-t", "--dirichlet-tol", help="tolerance to the edge", type=float, default=2e-2)
  parser.add_argument("--center-radius", help="radius of the 'center' dirichlet circle, relative to max xy extent", type=float, default=5e-2)
  parser.add_argument("--merge-tol", help="merge nodes tolerance", type=float, default=1e-6)
  parser.add_argument("--dirichlet", help="borders to clamp as dirichlet",
                      nargs="+", default=["xmin=0b111111","xmax=0b111111"])
  parser.add_argument("--min-comp-size", type=int, default=10)
  parser.add_argument("--grid", type=int, nargs="+", metavar="N",
    help="generate grid graph, 1 arg: NxN, 2 args: NXxNY")
  parser.add_argument("--hex", type=int, nargs="+", metavar="N",
    help="generate hexagonal honeycomb graph, 1 arg: NxN, 2 args: NXxNY")
  parser.add_argument("--mikado", type=float, nargs="+", metavar="X",
    help="generate random mikado fiber graph on the unit square, "
         "1 arg: MASS (total fiber length), 2 args: MASS R (fiber length, default 0.05)")
  parser.add_argument("--seed", type=int, default=0,
    help="random seed for --mikado")
  parser.add_argument("--min-edge", type=float, default=None, metavar="LEN",
    help="minimum edge length for --mikado: nodes closer than LEN are merged "
         "(default: r*1e-4)")
  parser.add_argument("--clamp-xy", type=float, nargs="+", metavar="X", default=None,
    help="clamp network to xy bounding box, drop edges with any endpoint outside, relative, at most two args")
  parser.add_argument("--no-props", action="store_true",
                    help="do not write per-edge properties to h5")
  parser.add_argument("--rescale-bbox", action="store_true",
                    help="rescale network so xy bbox is 1x1 (z scaled by same factor)")
  parser.add_argument("--rescale-props", default=None,
                    help="rescale network material properties, format '1,2,3,...'")
  parser.add_argument("--quirk", default=None, choices=Network.QUIRKS,
                    help="apply quirk")
  parser.add_argument("--prop-cutoff", type=float, metavar="PCT", default=None,
                    help="floor each stiffness component (edgeProps cols 1..6) at its PCT-th "
                         "percentile to shrink coefficient contrast, raising at most ~PCT%% of "
                         "fibers per component. Prints a per-decade log histogram, reference "
                         "thresholds, and the applied floor, then writes the clipped network. "
                         "Applied to the final (post-clamp) fibers.")
  args = parser.parse_args()

  if args.rescale_props is not None:
    rescale_props = np.array(list(map(float, args.rescale_props.split(","))))

  network = Network()
  if args.grid is not None:
    if len(args.grid) == 1:
      nx = ny = args.grid[0]
    elif len(args.grid) == 2:
      nx, ny = args.grid
    else:
      parser.error("--grid takes 1 or 2 arguments")
    network.generate_grid(nx, ny)
    network.generate_synthetic_properties(width=0.1 / max(nx, ny))
  elif args.hex is not None:
    if len(args.hex) == 1:
      nx = ny = args.hex[0]
      nz = 1
    elif len(args.hex) == 2:
      nx, ny = args.hex
      nz = 1
    elif len(args.hex) == 3:
      nx, ny, nz = args.hex
    else:
      parser.error("--hex takes 1 or 2 arguments")
    network.generate_honeycomb(nx, ny, nz)
    network.generate_synthetic_properties(width=0.1 / max(nx, ny))
  elif args.mikado is not None:
    if len(args.mikado) == 1:
      mass, r = args.mikado[0], 0.05
    elif len(args.mikado) == 2:
      mass, r = args.mikado
    else:
      parser.error("--mikado takes 1 or 2 arguments")
    network.generate_mikado(mass, r=r, seed=args.seed, min_edge=args.min_edge)
    network.generate_synthetic_properties(width=1/mass)
  else:
    network.read_morgan(args.input, rescale_props=args.rescale_props, quirk=args.quirk)
    network.verify_nonzero()
    if args.clamp_xy is not None:
      if len(args.clamp_xy) == 1:
        fx = fy = args.clamp_xy[0]
      elif len(args.clamp_xy) == 2:
        fx, fy = args.clamp_xy
      else:
        parser.error("--clamp-xy takes 1 or 2 arguments")
      network.clamp_xy(fx, fy)
    network.node_edge_dedupe(args.merge_tol)
  tprint("info", network.info)
  if args.prop_cutoff is not None:
    network.prop_cutoff(args.prop_cutoff)
  if args.rescale_bbox:
    network.rescale_bbox()
  network.compute_types(args.dirichlet_tol)
  if args.grid is None and args.hex is None and args.mikado is None:
    network.drop_floating_and_small_components(args.min_comp_size)
  network.write_h5(args.output, no_props=args.no_props)
  network.write_vtkhdf_view(args.output)
