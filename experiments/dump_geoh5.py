#!/usr/bin/env python

import numpy as np
import time
import argparse
import pandas
import h5py
import yaml

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)

parser = argparse.ArgumentParser(description="dump_geoh5 by Joseph Holten")
parser.add_argument("-i", help="input")
args = parser.parse_args()

tprint("reading h5 file")
with h5py.File(args.i, "r") as f:
  g = f["domain"]
  types_points = g["types_points"][:]
  count_dir = types_points.sum()
  frac_dir = count_dir / len(types_points)
  tprint("count dir", count_dir)
  tprint(f"frac dir {frac_dir:.5e}")
