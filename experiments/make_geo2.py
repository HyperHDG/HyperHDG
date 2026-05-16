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


  def read_morgan(self, path):
    tprint("reading nodes")
    nodes   = pandas.read_csv(path + "/nodes.csv")
    nodes   = nodes.to_numpy()[:,1:]

    tprint("reading edges")
    edges   = pandas.read_csv(path + '/edges.csv')
    edges   = edges.to_numpy()[:,1:]


    tprint("reading edgeProps")
    edgeProps   = pandas.read_csv(path + '/edgeProperties.csv')
    edgeProps   = edgeProps.to_numpy()[:,2:]

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
        for unit, value in zip(
            ["length", "time", "weight"],
            f.readlines()
        ):
          info["unit_"+unit] = value.strip()
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
    assert(edgeProps_dim == 12)

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


  def write_h5(self, out):
    tprint(f"writing h5 file to '{out}'")
    with h5py.File(out + ".geo.h5", "w") as f:
      g = f.create_group("domain")
      g.create_dataset("points", data=self.nodes, compression="gzip")
      g.create_dataset("edges", data=self.edges, compression="gzip")
      if hasattr(self, "edgeProps") and self.edgeProps is not None:
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
    path = out + ".geo.h5"
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
  parser.add_argument("--clamp-xy", type=float, nargs="+", metavar="X", default=None,
    help="clamp network to xy bounding box, drop edges with any endpoint outside, relative, at most two args")
  args = parser.parse_args()

  network = Network()
  if args.grid is not None:
    if len(args.grid) == 1:
      nx = ny = args.grid[0]
    elif len(args.grid) == 2:
      nx, ny = args.grid
    else:
      parser.error("--grid takes 1 or 2 arguments")
    network.generate_grid(nx, ny)
  else:
    network.read_morgan(args.input)
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
  network.compute_types(args.dirichlet_tol)
  if args.grid is None:
    network.drop_floating_and_small_components(args.min_comp_size)
  network.write_h5(args.output)
  network.write_vtkhdf_view(args.output)
