#!/usr/bin/env python
import h5py
import numpy as np
import argparse

p = argparse.ArgumentParser()
p.add_argument("input", help="path to .vtkhdf file")
p.add_argument("-c", "--component", help="component")

args = p.parse_args()

comps = {'n': 0, "m": 3, "u": 6, "r": 9, "v": 12, "s": 15}

with h5py.File(args.input, "r") as f:
  values = f["VTKHDF/PointData/values"][:]
  npts = f["VTKHDF/NumberOfPoints"][0]
  nsteps = values.shape[0] // npts
  steps = values.reshape(nsteps, npts, -1)
  comps_start = 0
  comps_end = steps[0].shape[1]

if args.component:
  comp = comps[args.component]
  comps_start = comp
  comps_end = comp+3
else:
  comps_start = 0
  comps_end = steps[0].shape[1]

print("k,c,max,mean,ratio")
for k in range(nsteps):
    x = steps[k]
    for c in range(comps_start, comps_end):
        m = x[:, c].max()
        a = x[:, c].mean()
        r = m / a if a > 0 else 0.0
        print(f"{k:3d},{c:2d},{m:.3e},{a:.3e},{r:.3e}")
