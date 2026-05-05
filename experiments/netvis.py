#!/usr/bin/env pvpython
import argparse
from paraview.simple import *
from matplotlib.colors import to_rgb
from lxml import etree
import h5py
import tempfile
import os
import sys
import copy

VIEWS = {
  "top":    {"position": (0, 0, 1), "focal": (0, 0, 0), "up": (0, 1, 0), "parallel": 1.0},
  "bottom": {"position": (0, 0,-1), "focal": (0, 0, 0), "up": (1, 0, 0), "parallel": 1.0},
  "side":   {"position": (0, -1, 0), "focal": (0, 0, 0), "up": (0, 0, 1), "parallel": 1.0},
  "iso":   {"position": (1, 1, 1), "focal": (0, 0, 0), "up": (0, 0, 1), "parallel": 1.0},
}

def apply_view(cam, v):
  cam.SetPosition(*v["position"])
  cam.SetFocalPoint(*v["focal"])
  cam.SetViewUp(*v["up"])
  cam.SetParallelScale(v["parallel"])

def write_xdmf3(path, domain, partition=None, solution=None, h5_i64=False, trace=None):
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
        tp_size = f["domain/types_points"].dtype.itemsize

    static_children = []

    topo = etree.Element("Topology", TopologyType="Polyline", NodesPerElement="2",
                         NumberOfElements=str(n_edges))
    etree.SubElement(topo, "DataItem", Format="HDF", DataType="Int",
                     Dimensions=f"{n_edges} 2").text = f"{domain}:/domain/edges"
    static_children.append(topo)

    geo = etree.Element("Geometry", GeometryType="XYZ")
    etree.SubElement(geo, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                     Dimensions=f"{n_points} 3").text = f"{domain}:/domain/points"
    static_children.append(geo)

    attr = etree.Element("Attribute", Name="properties", Center="Cell",
                         AttributeType="Matrix")
    etree.SubElement(attr, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                     Dimensions=f"{n_edges} {n_props}").text = f"{domain}:/domain/properties"
    static_children.append(attr)

    attr = etree.Element("Attribute", Name="types_points", Center="Node")
    etree.SubElement(attr, "DataItem", Format="HDF", DataType="Int", Precision=str(tp_size),
                     Dimensions=str(n_points)).text = f"{domain}:/domain/types_points"
    static_children.append(attr)

    if partition:
        partition = os.path.abspath(partition)
        attr = etree.Element("Attribute", Name="partition", Center="Node")
        etree.SubElement(attr, "DataItem", Format="HDF", DataType="Int", Precision=str(isize),
                         Dimensions=str(n_points)).text = f"{partition}:/net2as_part"
        static_children.append(attr)

    if solution:
        solution = os.path.abspath(solution)
        with h5py.File(solution, "r") as f:
            dset = list(f.keys())[0]
            srows, scols = f[dset].shape
            size = f[dset].dtype.itemsize
            assert srows == n_points, f"/{dset}: expected n_points == {n_points} rows, found {srows}"
            assert scols == 6, f"/{dset}: expected 6 cols, found {scols}"
            assert size == fsize, f"/{dset}: expected size {fsize} data type, found {size}"
        attr = etree.Element("Attribute", Name="solution", Center="Node",
                             AttributeType="Matrix")
        etree.SubElement(attr, "DataItem", Format="HDF", DataType="Float", Precision=str(fsize),
                         Dimensions=f"{srows} {scols}").text = f"{solution}:/{dset}"
        static_children.append(attr)

    xdmf = etree.Element("Xdmf", Version="3.0")
    dom = etree.SubElement(xdmf, "Domain")

    if trace:
        trace = os.path.abspath(trace)
        grid = etree.SubElement(dom, "Grid", Name="network",
                                GridType="Collection", CollectionType="Temporal")
        with h5py.File(trace, "r") as f:
            group = f["trace"]
            timesteps = sorted(k for k in group.keys() if k.startswith("timestep_"))
            for k, ts in enumerate(timesteps):
                d = group[ts]
                trows, tcols = d.shape
                assert d.dtype.itemsize == fsize
                sub = etree.SubElement(grid, "Grid", Name=ts, GridType="Uniform")
                etree.SubElement(sub, "Time", Value=str(k))
                for child in static_children:
                    sub.append(copy.deepcopy(child))
                attr = etree.SubElement(sub, "Attribute", Name="trace",
                                        Center="Node", AttributeType="Matrix")
                etree.SubElement(attr, "DataItem", Format="HDF", DataType="Float",
                                 Precision=str(fsize),
                                 Dimensions=f"{trows} {tcols}"
                                 ).text = f"{trace}:/trace/{ts}"
    else:
        grid = etree.SubElement(dom, "Grid", Name="network", GridType="Uniform")
        for child in static_children:
            grid.append(child)

    etree.ElementTree(xdmf).write(path, xml_declaration=True, pretty_print=True)
    return xdmf

def netvis(domain, partition=None, radius=None, use_tubes=True, output=None, show=True, resolution=(1000, 1000), solution=None, dirichlet=False, no_ref=False, trace=None, duration=5):
  if radius is None:
    try:
      with h5py.File(domain, "r") as f:
        radius = f["domain"].attrs["radius"]
    except KeyError as e:
      print("ERROR: no radius provided, and couldnt find radius in the domain file", file=sys.stderr)
  else:
    with h5py.File(domain, "a") as f:
      f["domain"].attrs["radius"] = radius

  print(f"xmf3_path: {args.xmf3}")
  write_xdmf3(args.xmf3, domain, partition=partition, solution=solution, trace=trace)

  pipe = Xdmf3ReaderS(FileName=[args.xmf3])
  pipe.UpdatePipeline()

  if trace:
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.PlayMode = "Snap To TimeSteps"

  pipe = ExtractSurface(Input=pipe)
  pipe.UpdatePipeline()

  disp_src = "solution" if solution else ("trace" if trace else None)

  if disp_src:
    pipe = Calculator(Input=pipe)
    pipe.AttributeType = "Point Data"
    pipe.ResultArrayName = "displacement"
    pipe.Function = f"{disp_src}_0*iHat + {disp_src}_1*jHat + {disp_src}_2*kHat"
    pipe.UpdatePipeline()

    if not no_ref:
      outline = Outline(Input=pipe)
      ref = Show(outline, GetActiveView())
      ref.AmbientColor = list(to_rgb(args.fg))
      ref.DiffuseColor = list(to_rgb(args.fg))
      ref.Opacity = 0.4

    pipe = WarpByVector(Input=pipe)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = 1.0
    pipe.UpdatePipeline()

  if use_tubes:
    pipe = Tube(Input=pipe)
    pipe.Radius = radius
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
  cam = GetActiveCamera()
  apply_view(cam, VIEWS[args.view])
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
  p.add_argument("-t", "--trace", help="trace data .sol.h5", default=None)
  p.add_argument("--fg", default="white", help="fg color")
  p.add_argument("--bg", default="black", help="bg color")
  p.add_argument("--xmf3", default=tempfile.mktemp(suffix=".xmf3"), help="xmf3 path")
  p.add_argument("--no-tubes", default=False, help="select tubes", action="store_true")
  p.add_argument("--h5-i64", default=False, help="select 64-bit integers in h5", action="store_true")
  p.add_argument("-o", "--output", default="netvis.png", help="save screenshot of network (after interact)")
  p.add_argument("--resolution", default="1000x1000", help="render resolution")
  p.add_argument("--duration", type=float, default=5, help="animation duration")
  p.add_argument("--dirichlet", default=False, help="color Dirichlet nodes (types_points == 1)", action="store_true")
  p.add_argument("--view", choices=list(VIEWS) + [None], default="top",
                 help="named camera view")
  p.add_argument("--no-ref", help="no reference view", action='store_true')
  args = p.parse_args()

  resolution = map(int, args.resolution.split("x"))
  netvis(args.domain, partition=args.partition, radius=args.radius, resolution=resolution, output=args.output, solution=args.solution, dirichlet=args.dirichlet, no_ref=args.no_ref, trace=args.trace, duration=args.duration)
