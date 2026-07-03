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
    tprint(f"generating honeycomb graph {nx} x {ny} x {nz}")
    s3 = np.sqrt(3.0)
    basis = np.array([
        [0.0,  0.0     ],
        [1.0,  0.0     ],
        [1.5,  0.5*s3  ],
        [2.5,  0.5*s3  ],
    ])
    a1 = np.array([3.0, 0.0])
    a2 = np.array([0.0, s3 ])
    ii, jj = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    origins = jj[..., None] * a1 + ii[..., None] * a2
    pts = origins[..., None, :] + basis[None, None, :, :]
    nodes2 = pts.reshape(-1, 2)
    def idx(i, j, k): return (i * nx + j) * 4 + k
    e = []
    i, j = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    e.append(np.stack([idx(i,j,0).ravel(), idx(i,j,1).ravel()], 1))
    e.append(np.stack([idx(i,j,1).ravel(), idx(i,j,2).ravel()], 1))
    e.append(np.stack([idx(i,j,2).ravel(), idx(i,j,3).ravel()], 1))
    i, j = np.meshgrid(np.arange(ny), np.arange(nx-1), indexing='ij')
    e.append(np.stack([idx(i,j,3).ravel(), idx(i,j+1,0).ravel()], 1))
    i, j = np.meshgrid(np.arange(ny-1), np.arange(nx), indexing='ij')
    e.append(np.stack([idx(i,j,2).ravel(), idx(i+1,j,1).ravel()], 1))
    i, j = np.meshgrid(np.arange(1, ny), np.arange(1, nx), indexing='ij')
    e.append(np.stack([idx(i,j,0).ravel(), idx(i-1,j-1,3).ravel()], 1))
    edges2d = np.vstack(e).astype(np.int64)

    # Drop left-boundary atom-0 and right-boundary atom-3 stubs
    n2d = nodes2.shape[0]
    keep = np.ones(n2d, dtype=bool)
    i_all = np.arange(ny)
    keep[(i_all * nx + 0) * 4 + 0] = False
    keep[(i_all * nx + (nx-1)) * 4 + 3] = False
    edge_keep = keep[edges2d[:,0]] & keep[edges2d[:,1]]
    edges2d = edges2d[edge_keep]
    remap = np.full(n2d, -1, dtype=edges2d.dtype)
    remap[keep] = np.arange(keep.sum())
    edges2d = remap[edges2d]
    nodes2 = nodes2[keep]

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

    sides = {
      "xmin": nodes[:,0] - mins[0] <= tol * dims[0],
      "xmax": maxs[0] - nodes[:,0] <= tol * dims[0],
      "ymin": nodes[:,1] - mins[1] <= tol * dims[1],
      "ymax": maxs[1] - nodes[:,1] <= tol * dims[1]
    }

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
  parser.add_argument("--merge-tol", help="merge nodes tolerance", type=float, default=1e-6)
  parser.add_argument("--dirichlet", help="borders to clamp as dirichlet",
                      nargs="+", default=["xmin=0b111111","xmax=0b111111"])
  parser.add_argument("--min-comp-size", type=int, default=10)
  parser.add_argument("--grid", type=int, nargs="+", metavar="N",
    help="generate grid graph, 1 arg: NxN, 2 args: NXxNY")
  parser.add_argument("--hex", type=int, nargs="+", metavar="N",
    help="generate hexagonal honeycomb graph, 1 arg: NxN, 2 args: NXxNY")
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
  if args.rescale_bbox:
    network.rescale_bbox()
  network.compute_types(args.dirichlet_tol)
  if args.grid is None and args.hex is None:
    network.drop_floating_and_small_components(args.min_comp_size)
  network.write_h5(args.output, no_props=args.no_props)
  network.write_vtkhdf_view(args.output)
