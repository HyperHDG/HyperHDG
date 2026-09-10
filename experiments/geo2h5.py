#!/usr/bin/env python

"""
geo2h5.py -- convert a legacy ASCII .geo domain into a .geo.h5 file.

The .geo format (see domains/*.geo, parsed by read_domain_geo() in
include/HyperHDG/read_domain.hxx) and the HDF5 format (read_domain_hdf5(),
written by make_geo2.py) carry the same information, but only the latter is
understood by the python tooling (netvis.py, gortz_constants.py, ...) and by
the distributed reader.  This script bridges the two.

The HDF5 layout mirrors make_geo2.py exactly:

    /domain/points        (n_points, space_dim)  f8
    /domain/edges         (n_edges,  2^hyEdge_dim) i4   -- point indices
    /domain/types_faces   (n_edges,  2*hyEdge_dim) i4
    /domain/types_points  (n_points, 1)          i4     -- derived, see below
    /domain/properties    (n_edges,  n_props)    f8     -- if present
    /VTKHDF/...                                         -- view for ParaView

Two restrictions of the HDF5 reader are checked here rather than silently
mangled: read_domain_hdf5() fills both hyNodes_hyEdge and points_hyEdge from
the single 'edges' dataset and allocates n_hyNodes == n_points, so the .geo
file must number its hypernodes exactly like its points.  For the
hyEdge_dim == 1 network files this always holds.

types_points has no .geo counterpart; it is derived as the maximum face type
seen at each hypernode, which reproduces the Dirichlet marker make_geo2.py
writes (a node is Dirichlet as soon as one incident face is).

Examples
--------
    geo2h5.py domains/cross.geo                  # -> domains/cross.geo.h5
    geo2h5.py domains/line3.geo -o /tmp/l3.geo.h5

Author
------
Joseph Holten, KIT, 2026.
"""

import argparse
import os
import sys

import h5py
import numpy as np


def clean(line):
  """Strip a '#' comment and the trailing ';' the header lines carry."""
  return line.split("#")[0].replace(";", " ").strip()


class GeoFile:
  """Sequential reader over the non-empty, non-comment lines of a .geo file."""

  def __init__(self, path):
    with open(path) as f:
      self.lines = [c for c in map(clean, f) if c]
    self.path = path
    self.pos = 0

  def seek(self, keyword):
    """Advance to the line starting with `keyword`; return its remainder."""
    while self.pos < len(self.lines):
      line = self.lines[self.pos]
      self.pos += 1
      if line.split()[0] == keyword:
        return line.split(None, 1)[1] if len(line.split()) > 1 else ""
    raise ValueError(f"keyword '{keyword}' not found in {self.path}")

  def scalar(self, keyword):
    """Read a 'Keyword = value;' header entry."""
    rest = self.seek(keyword).split()
    if not rest or rest[0] != "=":
      raise ValueError(f"'{keyword}' is not followed by '=' in {self.path}")
    return int(rest[1])

  def block(self, keyword, rows, cols, dtype):
    """Seek `keyword`, then read `rows` lines of `cols` numbers each."""
    self.seek(keyword)
    return self.rows(keyword, rows, cols, dtype)

  def rows(self, keyword, rows, cols, dtype):
    """Read `rows` lines of `cols` numbers each from the current position."""
    out = np.zeros((rows, cols), dtype=dtype)
    for i in range(rows):
      if self.pos >= len(self.lines):
        raise ValueError(f"{self.path}: '{keyword}' block ends after {i}/{rows} rows")
      vals = self.lines[self.pos].split()
      self.pos += 1
      if len(vals) > cols:
        raise ValueError(f"{self.path}: '{keyword}' row {i} has {len(vals)} > {cols} entries")
      # short rows are zero-padded, matching the C++ reader: a failed
      # istream extraction leaves the target at 0 (e.g. cross.geo declares
      # Space_Dim = 3 but lists only two coordinates per point)
      out[i, :len(vals)] = [dtype(v) for v in vals]
    return out


def read_geo(path):
  geo = GeoFile(path)

  space_dim = geo.scalar("Space_Dim")
  hyEdge_dim = geo.scalar("HyperEdge_Dim")
  n_points = geo.scalar("N_Points")
  n_hyNodes = geo.scalar("N_HyperNodes")
  n_hyEdges = geo.scalar("N_HyperEdges")

  n_faces = 2 * hyEdge_dim
  n_verts = 1 << hyEdge_dim

  points = geo.block("POINTS:", n_points, space_dim, float)
  hyNodes = geo.block("HYPERNODES_OF_HYPEREDGES:", n_hyEdges, n_faces, int)
  types_faces = geo.block("TYPES_OF_HYPERFACES:", n_hyEdges, n_faces, int)
  edges = geo.block("POINTS_OF_HYPEREDGES:", n_hyEdges, n_verts, int)

  props = None
  try:
    head = geo.seek("HYPEREDGE_PROPERTIES:")
  except ValueError:
    pass
  else:
    # the keyword line carries the property count: "HYPEREDGE_PROPERTIES: 12"
    props = geo.rows("HYPEREDGE_PROPERTIES:", n_hyEdges, int(head.split()[0]), float)

  if n_hyNodes != n_points:
    sys.exit(f"error: {path} has {n_hyNodes} hypernodes but {n_points} points; the h5 format "
             "assumes one hypernode per point")
  if not np.array_equal(hyNodes, edges):
    sys.exit(f"error: {path} numbers hypernodes differently from points; the h5 format stores "
             "a single 'edges' dataset for both")

  return dict(space_dim=space_dim, hyEdge_dim=hyEdge_dim, points=points,
              edges=edges, types_faces=types_faces, properties=props)


def types_points(n_points, edges, types_faces):
  """Per-node type = strongest (max) type among the faces meeting at that node."""
  out = np.zeros((n_points, 1), dtype=np.int32)
  np.maximum.at(out[:, 0], edges.ravel(), types_faces.ravel())
  return out


def write_h5(out, geo, vtkhdf=True):
  points = np.ascontiguousarray(geo["points"], dtype=np.float64)
  edges = np.ascontiguousarray(geo["edges"], dtype=np.int32)
  faces = np.ascontiguousarray(geo["types_faces"], dtype=np.int32)
  gz = dict(compression="gzip", compression_opts=1)

  with h5py.File(out, "w") as f:
    g = f.create_group("domain")
    g.create_dataset("points", data=points, **gz)
    g.create_dataset("edges", data=edges, **gz)
    g.create_dataset("types_faces", data=faces, **gz)
    g.create_dataset("types_points", data=types_points(len(points), edges, faces), **gz)
    if geo["properties"] is not None:
      g.create_dataset("properties",
                       data=np.ascontiguousarray(geo["properties"], dtype=np.float64), **gz)
    g.attrs["size"] = points.max(axis=0) - points.min(axis=0)

    if vtkhdf:
      if geo["hyEdge_dim"] != 1:
        print("warning: VTKHDF view only supports hyEdge_dim == 1, skipping", file=sys.stderr)
        return
      n_cells = edges.shape[0]
      n_conn = 2 * n_cells

      root = f.create_group("VTKHDF")
      root.attrs.create("Version", [2, 0], dtype="int64")
      root.attrs.create("Type", np.bytes_("UnstructuredGrid"))
      root["Points"] = h5py.SoftLink("/domain/points")

      layout = h5py.VirtualLayout(shape=(n_conn,), dtype=edges.dtype)
      layout[...] = h5py.VirtualSource(out, "domain/edges", shape=edges.shape)
      root.create_virtual_dataset("Connectivity", layout)

      root.create_dataset("Offsets", data=np.arange(0, n_conn + 2, 2, dtype=np.int64), **gz)
      root.create_dataset("Types", data=np.full(n_cells, 3, dtype=np.uint8), **gz)
      root.create_dataset("NumberOfPoints", data=np.array([len(points)], dtype=np.int64))
      root.create_dataset("NumberOfCells", data=np.array([n_cells], dtype=np.int64))
      root.create_dataset("NumberOfConnectivityIds", data=np.array([n_conn], dtype=np.int64))

      pd = root.create_group("PointData")
      pd["types_points"] = h5py.SoftLink("/domain/types_points")
      if geo["properties"] is not None:
        cd = root.create_group("CellData")
        cd["properties"] = h5py.SoftLink("/domain/properties")


def main():
  p = argparse.ArgumentParser(description="convert a .geo domain to .geo.h5")
  p.add_argument("input", help="input .geo file")
  p.add_argument("-o", "--output", help="output .geo.h5 file (default: input with .geo.h5)")
  p.add_argument("--no-vtkhdf", action="store_true",
                 help="omit the VTKHDF view (netvis.py/ParaView then cannot read the file)")
  args = p.parse_args()

  out = args.output or os.path.splitext(args.input)[0] + ".geo.h5"
  geo = read_geo(args.input)
  write_h5(out, geo, vtkhdf=not args.no_vtkhdf)

  print(f"{args.input} -> {out}: {len(geo['points'])} points, {len(geo['edges'])} edges, "
        f"space_dim {geo['space_dim']}, hyEdge_dim {geo['hyEdge_dim']}, "
        f"{0 if geo['properties'] is None else geo['properties'].shape[1]} properties")


if __name__ == "__main__":
  main()
