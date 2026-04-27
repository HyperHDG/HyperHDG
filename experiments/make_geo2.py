#!/usr/bin/env python

import numpy as np
import time
import argparse
import pandas
import h5py
import yaml

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)


parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("-i", help="input", default=".")
parser.add_argument("-o", help="output", default="graph")
parser.add_argument("-t", help="tolerance to the edg", type=float, default=2e-2)
parser.add_argument("--dirichlet", help="dimension to clamp outer most as dirichlet", nargs="+", type=int, default=[0, 1])
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

