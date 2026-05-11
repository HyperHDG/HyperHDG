#!/usr/bin/env pvpython
import argparse
import sys
import os
from paraview.simple import *
from matplotlib.colors import to_rgb

class View:
  VIEWS = {
    "top":    {"position": (0, 0, 1), "focal": (0, 0, 0), "up": (0, 1, 0)},
    "bottom": {"position": (0, 0,-1), "focal": (0, 0, 0), "up": (1, 0, 0)},
    "side":   {"position": (0,-1, 0), "focal": (0, 0, 0), "up": (0, 0, 1)},
    "iso":    {"position": (1, 1, 1), "focal": (0, 0, 0), "up": (0, 0, 1)},
  }

  def __init__(self, view):
    self.view = View.VIEWS[view]

  def orient(self, cam):
    v = self.view
    cam.SetPosition(*v["position"])
    cam.SetFocalPoint(*v["focal"])
    cam.SetViewUp(*v["up"])

def has_point_array(proxy, name):
  pdi = proxy.GetPointDataInformation()
  for i in range(pdi.GetNumberOfArrays()):
    if pdi.GetArray(i).GetName() == name:
      return True
  return False


# --- Pipeline ops ------------------------------------------------------------
# Each op is a callable: pipe -> pipe. May skip and return input unchanged.

class Warp:
  def __init__(self, components=(6, 7, 8), scale=1.0, source_array="values"):
    self.components = components
    self.scale = scale
    self.source_array = source_array

  def __call__(self, pipe):
    if not has_point_array(pipe, self.source_array):
      print(f"warning: Warp: '{self.source_array}' not found, skipping",
            file=sys.stderr)
      return pipe
    cx, cy, cz = self.components
    calc = Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "displacement"
    calc.Function = (f"{self.source_array}_{cx}*iHat + "
                     f"{self.source_array}_{cy}*jHat + "
                     f"{self.source_array}_{cz}*kHat")
    calc.UpdatePipeline()
    pipe = WarpByVector(Input=calc)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = self.scale
    pipe.UpdatePipeline()
    return pipe


# --- Display config ----------------------------------------------------------
# Applied to the Show() proxy after the pipeline is rendered.

class SolidColor:
  def __init__(self, color="white"):
    self.color = color

  def __call__(self, display, pipe):
    display.SetScalarColoring(None, 0)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb


# --- Runner ------------------------------------------------------------------

def netvis(path, ops=(),
           fg=SolidColor("white"), bg="black", view="iso", resolution=(1000, 1000),
           axis=True, output=None, show=True):
  if not os.path.isfile(path):
    sys.exit(f"error: file not found: {path}")

  pipe = VTKHDFReader(FileName=[path])
  pipe.UpdatePipeline()
  for op in ops:
    pipe = op(pipe)

  rview = GetActiveViewOrCreate("RenderView")
  display = Show(pipe, rview)
  display.Representation = "Wireframe"
  fg(display, pipe)

  rview.ViewSize = list(resolution)
  rview.Background = list(to_rgb(bg))
  rview.UseColorPaletteForBackground = 0
  rview.OrientationAxesVisibility = int(axis)

  view.orient(GetActiveCamera())
  ResetCamera()
  Render()
  if show:
    Interact()
  if output:
    SaveScreenshot(output, rview, TransparentBackground=1)


# --- CLI ---------------------------------------------------------------------

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--fg", default="white")
  p.add_argument("--bg", default="black")
  p.add_argument("--view", choices=list(View.VIEWS), default="top")
  p.add_argument("--resolution", default="1000x1000")
  p.add_argument("--warp", action="store_true",
                 help="warp by displacement (components 6-8 of 'values')")
  p.add_argument("--warp-scale", type=float, default=1.0)
  p.add_argument("-o", "--output", default=None, help="optional screenshot path")
  p.add_argument("--no-show", action="store_true", help="skip Interact()")
  p.add_argument("--no-axis", action="store_true")
  args = p.parse_args()

  resolution = tuple(map(int, args.resolution.split("x")))

  ops = []
  if args.warp:
    ops.append(Warp(scale=args.warp_scale))

  fg = SolidColor(args.fg)

  netvis(args.input, ops=ops, fg=fg, bg=args.bg, view=View(args.view), axis=not args.no_axis,
         resolution=resolution, output=args.output, show=not args.no_show)
