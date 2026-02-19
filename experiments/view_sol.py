#!/usr/bin/env python

import numpy as np
import argparse
import paraview.simple as pv
import math

parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("--network", help="input", default="network.vtk")
parser.add_argument("--save", help="save to file", default="network.png")
parser.add_argument("--noshow", help="noshow", default=False, action="store_true")
parser.add_argument("--data", help="which data to display, one of: binary file, SIN, --")
args = parser.parse_args()

pipe = pv.OpenDataFile(args.network)
pipe.UpdatePipeline()

bounds = np.array(pipe.GetDataInformation().GetBounds()).reshape(3,2).T
dims = bounds[1] - bounds[0]
print(dims)

if args.data == "SIN":
  pipe = pv.Calculator(Input=pipe)
  pipe.ResultArrayName = "result"
  pipe.Function = f"sin(coordsX / {dims[0]} * 2 * {math.pi}) * sin(coordsY / {dims[1]} * 2 * {math.pi})"

  pipe = pv.WarpByScalar(Input=pipe)
  pipe.Scalars = ['POINTS', 'result']
  pipe.ScaleFactor = np.max(dims)/10

pipe = pv.Show(pipe)
view = pv.GetActiveViewOrCreate("RenderView")
pv.Render()

if args.save: pv.SaveScreenshot(args.save)
if not args.noshow: pv.Interact()
