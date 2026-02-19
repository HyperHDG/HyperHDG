#!/usr/bin/env python

import sys
import numpy as np
import argparse
import paraview.simple as pv
import math
from vtk.numpy_interface import dataset_adapter as dsa
import vtk

parser = argparse.ArgumentParser(description="make_geo2 by Joseph Holten")
parser.add_argument("--network", help="input", default="network.vtk")
parser.add_argument("--save", help="save to file", default="network.png")
parser.add_argument("--noshow", help="noshow", default=False, action="store_true")
parser.add_argument("--data", help="which data to display, one of: binary file, SIN, --")
parser.add_argument("--float_t", help="which binary float type to use, e.g. float64", default="float64")
parser.add_argument("--data_dim", help="number of components in data, default: 6, first three used as x,y,z displacement",
                    default=6, type=int)
args = parser.parse_args()

if args.data == "SIN":
  print("NOTE: APPLYING DEBUG SIN DATA")

  pipe = pv.OpenDataFile(args.network)
  pipe.UpdatePipeline()

  n_points = pipe.GetDataInformation().GetNumberOfPoints()
  print("n_points", n_points)

  bounds = np.array(pipe.GetDataInformation().GetBounds()).reshape(3,2).T
  dims = bounds[1] - bounds[0]

  pipe = pv.Calculator(Input=pipe)
  pipe.ResultArrayName = "result"
  pipe.Function = f"sin(coordsX / {dims[0]} * 2 * {math.pi}) * sin(coordsY / {dims[1]} * 2 * {math.pi})"

  pipe = pv.WarpByScalar(Input=pipe)
  pipe.Scalars = ['POINTS', 'result']
  pipe.ScaleFactor = np.max(dims)/10

else:
  if args.data == "--":
    print("NOTE: READING DATA FROM STDIN")
    data = np.frombuffer(sys.stdin.buffer.read(), dtype=args.float_t)
  else:
    data = np.fromfile(args.data, dtype=args.float_t)
  print(data.shape)

  data = data.reshape(-1, args.data_dim)[:,:3]

  reader = vtk.vtkUnstructuredGridReader()
  reader.SetFileName(args.network)
  reader.Update()

  n_points = reader.GetOutput().GetNumberOfPoints()
  print("n_points", n_points)

  wrapped = dsa.WrapDataObject(reader.GetOutput())
  wrapped.PointData.append(data, "data")
  pipe = pv.TrivialProducer()
  pipe.GetClientSideObject().SetOutput(reader.GetOutput())
  pipe.UpdatePipeline()

  pipe = pv.WarpByVector(Input=pipe)
  pipe.Vectors = ['POINTS', 'data']

pipe = pv.Show(pipe)
view = pv.GetActiveViewOrCreate("RenderView")
pv.Render()

if args.save: pv.SaveScreenshot(args.save)
if not args.noshow: pv.Interact()
