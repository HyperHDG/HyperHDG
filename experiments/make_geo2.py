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


parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("-i", help="input", default=".")
parser.add_argument("-o", help="output", default="graph")
parser.add_argument("-t", help="tolerance to the edg", type=float, default=2e-2)
parser.add_argument("--merge-tol", help="merge nodes tolerance", type=float, default=1e-6)
parser.add_argument("--dirichlet", help="dimension to clamp outer most as dirichlet", nargs="+", type=int, default=[0, 1])
parser.add_argument("--min-comp-size", type=int, default=10)
args = parser.parse_args()

tprint("reading nodes")
nodes   = pandas.read_csv(args.i + "/nodes.csv")
nodes   = nodes.to_numpy()[:,1:]
n_nodes = nodes.shape[0]

tprint("reading edges")
edges   = pandas.read_csv(args.i + '/edges.csv')
edges   = edges.to_numpy()[:,1:]
n_edges = edges.shape[0]

tprint("reading edgeProps")
edgeProps   = pandas.read_csv(args.i + '/edgeProperties.csv')
edgeProps   = edgeProps.to_numpy()[:,2:]
n_edgeProps = edgeProps.shape[0]
edgeProps_dim = edgeProps.shape[1]

tprint("nodes", nodes.shape)
tprint("edges", edges.shape)
tprint("edgeProps", edgeProps.shape)

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

close_pairs = tree.query_pairs(args.merge_tol, output_type='ndarray') # (n_pairs, 2), i < j
tprint(f"merging {len(close_pairs)} pairs closer than {args.merge_tol:.3e}")

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

try:
  with open(args.i + '/info.txt') as f:
    info = yaml.safe_load(f)
except FileNotFoundError:
  info = {}

try:
  with open(args.i + '/units.txt') as f:
    for unit, value in zip(
        ["length", "time", "weight"],
        f.readlines()
    ):
      info["unit_"+unit] = value.strip()
except FileNotFoundError:
  units = {}


mins = nodes.min(axis=0)
maxs = nodes.max(axis=0)
dims = maxs - mins

info["size"] = dims

tprint("info", info)


d = args.dirichlet
types_points = np.where(
  np.any(
    ((nodes[:, d] - mins[d]) < args.t * dims[d]) |
    ((maxs[d] - nodes[:, d]) < args.t * dims[d]),
    axis=1
  ),
  1, 0
).astype(np.int32)

count_dir = types_points.sum()
frac_dir = count_dir / len(types_points)
tprint("count dir", count_dir)
tprint(f"frac dir {frac_dir:.5e}")

types_faces = types_points[edges].astype(np.int32)

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

keep_comp = sizes >= args.min_comp_size  # boolean per component
# drop floating components (no Dirichlet node) regardless of size:
keep_comp &= ~free_mask
keep_node = keep_comp[labels]
tprint(f"keeping {keep_comp.sum()}/{n_comp} components, {keep_node.sum()}/{n_nodes} nodes")

# remap: old node index -> new node index, -1 for dropped
remap = np.full(n_nodes, -1, dtype=edges.dtype)
remap[keep_node] = np.arange(keep_node.sum())

nodes        = nodes[keep_node]
types_points = types_points[keep_node]
n_nodes      = len(nodes)

keep_edge = keep_node[edges[:,0]] & keep_node[edges[:,1]]
edges     = remap[edges[keep_edge]]
edgeProps = edgeProps[keep_edge]
types_faces  = types_faces[keep_edge]
n_edges   = len(edges)
tprint(f"after pruning: {n_nodes} nodes, {n_edges} edges")

tprint("writing h5 file")
with h5py.File(args.o + ".geo.h5", "w") as f:
    g = f.create_group("domain")
    g.create_dataset("points", data=nodes, compression="gzip")
    g.create_dataset("edges", data=edges, compression="gzip")
    g.create_dataset("properties", data=edgeProps, compression="gzip")
    g.create_dataset("types_points", data=types_points, compression="gzip")
    g.create_dataset("types_faces", data=types_faces, compression="gzip")
    for k, v in info.items():
      g.attrs[k] = v

