#!/usr/bin/env pvpython
import argparse
from paraview.simple import *
from matplotlib.colors import to_rgb
from lxml import etree
import h5py
import tempfile
import os
import sys

def write_xdmf3(path, domain, partition=None, solution=None, h5_i64=False):
  domain = os.path.abspath(domain)
  with h5py.File(domain, "r") as f:
    n_edges, ncols = f["domain/edges"].shape
    isize = f["domain/edges"].dtype.itemsize
    assert ncols == 2, f"domain/edges: expected 2 columns, found {ncols}"

    n_points, ncols = f["domain/points"].shape
    fsize = f["domain/points"].dtype.itemsize
    assert ncols == 3, f"domain/points: expected 3 columns, found {ncols}"

    nrows, n_props = f["domain/properties"].shape
    assert n_edges == nrows, \
      f"domain/properties: expected n_edges == {n_edges} rows, found {nrows}"

    tp_size = f["domain/types_points"].dtype.itemsize  # 4

  xdmf = etree.Element("Xdmf", Version="3.0")
  dom = etree.SubElement(xdmf, "Domain")
  grid = etree.SubElement(dom, "Grid", Name="network", GridType="Uniform")
  topo = etree.SubElement(grid, "Topology", TopologyType="Polyline", NodesPerElement="2",
                          NumberOfElements=str(n_edges))
  etree.SubElement(topo, "DataItem", Format="HDF", DataType="Int",
                   Dimensions=f"{n_edges} 2").text = f"{domain}:/domain/edges"
  geo = etree.SubElement(grid, "Geometry", GeometryType="XYZ")
  etree.SubElement(geo, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                   Dimensions=f"{n_points} 3").text = f"{domain}:/domain/points"

  attr = etree.SubElement(grid, "Attribute", Name="properties", Center="Cell",
                          AttributeType="Matrix")
  etree.SubElement(attr, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                   Dimensions=f"{n_edges} {n_props}").text = f"{domain}:/domain/properties"

  attr = etree.SubElement(grid, "Attribute", Name="types_points", Center="Node")
  etree.SubElement(attr, "DataItem", Format="HDF", DataType="Int", Precision=str(tp_size),
                   Dimensions=str(n_points)).text = f"{domain}:/domain/types_points"

  if partition:
    partition = os.path.abspath(partition)
    attr = etree.SubElement(grid, "Attribute", Name="partition", Center="Node")
    etree.SubElement(attr, "DataItem", Format="HDF", DataType="Int", Precision=str(isize),
                     Dimensions=str(n_points)).text = f"{partition}:/net2as_part"
  if solution:
    solution = os.path.abspath(solution)
    with h5py.File(solution, "r") as f:
      dset = list(f.keys())[0]
      nrows, ncols = f[dset].shape
      size = f[dset].dtype.itemsize
      assert nrows == n_points, f"/{dset}: expected n_points == {n_points} rows, found {nrows}"
      assert ncols == 6, f"/{dset}: expected 6 cols, found {ncols}"
      assert size == fsize, f"/{dset}: expected size {fsize} data type, found {size}"

    attr = etree.SubElement(grid, "Attribute", Name="solution", Center="Node",
                            AttributeType="Matrix")
    etree.SubElement(attr, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                     Dimensions=f"{nrows} {ncols}").text = f"{solution}:/{dset}"

  etree.ElementTree(xdmf).write(path, xml_declaration=True, pretty_print=True)

  return xdmf

def netvis(domain, partition=None, radius=None, use_tubes=True, output=None, show=True, resolution=(1000, 1000), solution=None, dirichlet=False):
  if radius is None:
    try:
      with h5py.File(domain, "r") as f:
        args.radius = f["domain"].attrs["radius"]
    except KeyError as e:
      print("ERROR: no radius provided, and couldnt find radius in the domain file", file=sys.stderr)
  else:
    with h5py.File(domain, "a") as f:
      f["domain"].attrs["radius"] = args.radius

  print(f"xmf3_path: {args.xmf3}")
  write_xdmf3(args.xmf3, domain, partition=partition, solution=solution)

  pipe = Xdmf3ReaderS(FileName=[args.xmf3])
  pipe.UpdatePipeline()

  pipe = ExtractSurface(Input=pipe)
  pipe.UpdatePipeline()

  if solution:
    pipe = Calculator(Input=pipe)
    pipe.AttributeType = "Point Data"
    pipe.ResultArrayName = "displacement"
    pipe.Function = "solution_0*iHat + solution_1*jHat + solution_2*kHat"
    pipe.UpdatePipeline()

    pipe = WarpByVector(Input=pipe)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = 1.0
    pipe.UpdatePipeline()

  if use_tubes:
    pipe = Tube(Input=pipe)
    pipe.Radius = args.radius
    pipe.NumberofSides = 4

  display = Show(pipe, GetActiveViewOrCreate("RenderView"))
  if partition:
    ColorBy(display, ("POINTS", "partition"))
    display.RescaleTransferFunctionToDataRange(True)
    lut = GetColorTransferFunction("partition")
    lut.ApplyPreset("Paired", True)
  elif dirichlet:
    ColorBy(display, ("POINTS", "types_points"))
    display.RescaleTransferFunctionToDataRange(True)
    lut = GetColorTransferFunction("types_points")
    lut.RGBPoints = [0.0, *to_rgb(args.fg), 1.0, *to_rgb("red")]
    lut.ColorSpace = "RGB"
    display.SetScalarBarVisibility(GetActiveView(), False)
  else:
    display.SetScalarColoring(None, 0)
    display.DiffuseColor = list(to_rgb(args.fg))

  view = GetActiveView()
  view.ViewSize = list(resolution)
  view.Background = list(to_rgb(args.bg))
  view.UseColorPaletteForBackground = 0
  view.ResetCamera()
  Render()
  Interact()
  if output:
    view.OrientationAxesVisibility = 0
    SaveScreenshot(output, view, TransparentBackground=1)

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("-d", "--domain", help="domain file .geo.h5", required=True)
  p.add_argument("-r", "--radius", default=None, help="radius value", type=float)
  p.add_argument("-p", "--partition", help="discrete data, e.g. partition .h5", default=None)
  p.add_argument("-s", "--solution", help="continuous data, e.g. solution .h5", default=None)
  p.add_argument("--fg", default="gray", help="fg color")
  p.add_argument("--bg", default="black", help="bg color")
  p.add_argument("--xmf3", default=tempfile.mktemp(suffix=".xmf3"), help="xmf3 path")
  p.add_argument("--no-tubes", default=False, help="select tubes", action="store_true")
  p.add_argument("--h5-i64", default=False, help="select 64-bit integers in h5", action="store_true")
  p.add_argument("-o", "--output", default="netvis.png", help="save screenshot of network (after interact)")
  p.add_argument("--resolution", default="1000x1000", help="render resolution")
  p.add_argument("--dirichlet", default=False, help="color Dirichlet nodes (types_points == 1)", action="store_true")
  args = p.parse_args()

  resolution = map(int, args.resolution.split("x"))
  netvis(args.domain, partition=args.partition, radius=args.radius, resolution=resolution, output=args.output, solution=args.solution, dirichlet=args.dirichlet)
