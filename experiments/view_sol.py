#!/usr/bin/env python

import numpy as np
import time
import argparse
import pandas
import h5py
import yaml
import paraview.simple as pv

parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("-i", help="input", default="network.vtk")
parser.add_argument("-o", help="output", default="graph.png")
parser.add_argument("-n", help="noshow", default=True, action="store_true")
args = parser.parse_args()

reader = pv.OpenDataFile(args.i)
display = pv.Show(reader)
view = pv.GetActiveViewOrCreate("RenderView")
pv.Render()
if args.o: pv.SaveScreenshot(args.o)
if not args.n: pv.Interact()
