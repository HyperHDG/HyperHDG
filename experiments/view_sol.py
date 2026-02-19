#!/usr/bin/env python

import numpy as np
import argparse
import paraview.simple as pv
import math

parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("-i", help="input", default="network.vtk")
parser.add_argument("-o", help="output", default="graph.png")
parser.add_argument("-n", help="noshow", default=False, action="store_true")
parser.add_argument("--sin", help="debug sin", default=True, action="store_true")
args = parser.parse_args()

reader = pv.OpenDataFile(args.i)
reader.UpdatePipeline()

bounds = np.array(reader.GetDataInformation().GetBounds()).reshape(3,2).T
dims = bounds[1] - bounds[0]
print(dims)

calc = pv.Calculator(Input=reader)
calc.ResultArrayName = "result"
calc.Function = f"sin(coordsX / {dims[0]} * 2 * {math.pi}) * sin(coordsY / {dims[1]} * 2 * {math.pi})"

display = pv.Show(calc)
view = pv.GetActiveViewOrCreate("RenderView")
pv.Render()

if args.o: pv.SaveScreenshot(args.o)
if not args.n: pv.Interact()
