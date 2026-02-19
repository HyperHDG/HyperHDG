#!/usr/bin/env python

import numpy as np
import time
import argparse
import pandas
import h5py
import yaml

parser = argparse.ArgumentParser(description="dump_geoh5 by Joseph Holten")
parser.add_argument("input", help="input")
args = parser.parse_args()

print("path:", args.input)

with h5py.File(args.input, "r") as f:
  g = f["domain"]
  print("shapes:")
  for key in ["points", "edges", "types_points", "types_faces"]:
    print(f"  {key}: {list(g[key].shape)}")
  types_points = g["types_points"][:]
  print("dir:")
  count_dir = types_points.sum()
  frac_dir = count_dir / len(types_points)
  print("  count:", count_dir)
  print(f"  frac: {frac_dir:.5e}")
  print("info:")
  for attr, val in g.attrs.items():
    print(f"  {attr}: {val}")
