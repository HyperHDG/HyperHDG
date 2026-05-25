#!/usr/bin/env pvpython
import argparse
import sys
import os
import paraview.simple as pv
from matplotlib.colors import to_rgb
from colorsys import hsv_to_rgb
import math
import numpy as np
import matplotlib.pyplot as plt

class View:
  VIEWS = {
    "top":    {"position": (0, 0, 1), "focal": (0, 0, 0), "up": (0, 1, 0)},
    "bottom": {"position": (0, 0,-1), "focal": (0, 0, 0), "up": (1, 0, 0)},
    "pside":  {"position": (0,-1, 0), "focal": (0, 0, 0), "up": (0, 0, 1)},
    "iside":  {"position": (-0.1128, 0.9824, 0.1486), "focal": (0, 0, 0), "up": (0, -0.1496, 0.9887)},
    "iso":    {"position": (1, 1, 1), "focal": (0, 0, 0), "up": (0, 0, 1)},
  }

  def __init__(self, view):
    self.view = view

  def orient(self, rview):
    cam = pv.GetActiveCamera()
    v = View.VIEWS[self.view]
    if self.view != "iso":
      print("using parallel projection")
      rview.CameraParallelProjection = 1
    cam.SetPosition(*v["position"])
    cam.SetFocalPoint(*v["focal"])
    cam.SetViewUp(*v["up"])

def find_array(pipe, name):
  for assoc, getter in (("POINTS", pipe.GetPointDataInformation),
                        ("CELLS",  pipe.GetCellDataInformation)):
    info = getter()
    if any(info.GetArray(i).GetName() == name for i in range(info.GetNumberOfArrays())):
      return assoc
  return None

# --- Pipeline ops ------------------------------------------------------------
# Each op is a callable: pipe -> pipe. May skip and return input unchanged.

class Warp:
  def __init__(self, components=(6, 7, 8), scale=1.0, source_array="values", normal=None):
    self.components = components
    self.scale = scale
    self.source_array = source_array
    self.normal = normal

  def apply(self, pipe, view):
    if find_array(pipe, self.source_array) != "POINTS":
      print(f"warning: Warp: '{self.source_array}' not found, skipping",
            file=sys.stderr)
      return pipe
    calc = pv.Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "displacement"
    if self.normal is None:
      cx, cy, cz = self.components
      calc.Function = (f"{self.source_array}_{cx}*iHat + "
                       f"{self.source_array}_{cy}*jHat + "
                       f"{self.source_array}_{cz}*kHat")
    else:
      cx, = self.components
      calc.Function = (f"{self.source_array}*{self.normal}")
    calc.UpdatePipeline()
    pipe = pv.WarpByVector(Input=calc)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = self.scale
    pipe.UpdatePipeline()
    return pipe

class Tubes:
  def __init__(self, radius=None, sides=4):
    self.radius = radius
    self.sides = sides

  def apply(self, pipe, view):
    pipe = pv.ExtractSurface(Input=pipe)
    pipe.UpdatePipeline()
    tube = pv.Tube(Input=pipe)
    if self.radius is not None:
      tube.Radius = self.radius
    tube.NumberofSides = self.sides
    tube.UpdatePipeline()
    return tube

class Reference:
  def __init__(self, color="white", opacity=1.0):
    self.color = color
    self.opacity = opacity

  def apply(self, pipe, view):
    outline = pv.Outline(Input=pipe)
    outline.UpdatePipeline()
    display = pv.Show(outline, view)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb
    display.Opacity = self.opacity
    return pipe

# --- Display config ----------------------------------------------------------
# Applied to the Show() proxy after the pipeline is rendered.

class SolidColor:
  def __init__(self, color="white"):
    self.color = color

  def apply(self, pipe, rview):
    display = pv.Show(pipe, rview)
    display.SetScalarColoring(None, 0)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb


class ArrayColor:
  def __init__(self, spec, fg="white", invert=False, categories=""):
    """'name' or 'name:N' -> (name, component_or_None)."""
    self.invert = invert
    self.fg = fg

    self.categories = []
    items = categories.split(",")
    for item in items:
      if ":" in item:
        a,b = item.split(":")
        self.categories.extend(range(a,b+1))
      elif len(item) > 0:
        self.categories.append(item)

    if ":" in spec:
      self.name, self.comp = spec.rsplit(":", 1)
      self.comp = int(self.comp)
    else:
      self.name = spec
      self.comp = None

  def apply(self, pipe, rview):
    display = pv.Show(pipe, rview)
    assoc = find_array(pipe, self.name)
    if assoc is None:
      print(f"warning: '{self.name}' not found", file=sys.stderr)
      return
    if self.comp is None:
      target = (assoc, self.name)
    else:
      target = (assoc, self.name, self.comp)
    pv.ColorBy(display, target)
    ctf = pv.GetColorTransferFunction(self.name)
    comp = self.comp if self.comp is not None else 0
    rng = display.GetArrayInformationForColorArray().GetComponentRange(comp)
    ctf.ApplyPreset("Cool to Warm", True)
    if len(self.categories) != 0:
      ctf.InterpretValuesAsCategories = 1
      ctf.AnnotationsInitialized = 1

      cats = [0] + [c for c in self.categories if c != 0]

      annotations = []
      for v in cats:
          label = f"{int(v):06b}"
          annotations.extend([str(v), label])
      ctf.Annotations = annotations

      colors = list(to_rgb(self.fg))
      n = len(cats) - 1
      for i in range(n):
          colors.extend(hsv_to_rgb(i / n, 0.7, 0.9))
      ctf.IndexedColors = colors
      ctf.IndexedOpacities = [1.0] * len(cats)
    else:
      ctf.ApplyPreset("Cool to Warm", True)
      M = max(abs(rng[0]), abs(rng[1]))
      ctf.RescaleTransferFunction(-M, M)
      if self.invert:
        ctf.InvertTransferFunction()
    display.SetScalarBarVisibility(pv.GetActiveView(), True)

class CoarseArrows:
  def __init__(self, disp_components=(6, 7, 8), rot_components=(9, 10, 11),
    source_array="values", resolution=(10, 10), plane="xy", kernel_radius=None,
    warp_scale=1.0, offset_z=0.0, scale=3.0, color="red"):
    self.disp_components = disp_components
    self.rot_components = rot_components
    self.source_array = source_array
    self.resolution = resolution       # 2D: (n1, n2) in-plane
    self.plane = plane                 # "xy", "xz", "yz"
    self.kernel_radius = kernel_radius # None -> auto from grid spacing
    self.warp_scale = warp_scale     # match the main Warp scale
    self.offset_z = offset_z         # fixed lift after warping
    self.scale = scale
    self.color = color

  def apply(self, pipe, rview):
    if find_array(pipe, self.source_array) != "POINTS":
      print(f"warning: CoarseArrows: '{self.source_array}' not found, skipping",
        file=sys.stderr)
      return pipe

    # build both vector fields on the (reference) input
    dx_, dy_, dz_ = self.disp_components
    rx, ry, rz = self.rot_components
    calc = pv.Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "rotation"
    calc.Function = (f"{self.source_array}_{rx}*iHat + "
                     f"{self.source_array}_{ry}*jHat + "
                     f"{self.source_array}_{rz}*kHat")
    calc.UpdatePipeline()

    calc2 = pv.Calculator(Input=calc)
    calc2.AttributeType = "Point Data"
    calc2.ResultArrayName = "displacement"
    calc2.Function = (f"{self.source_array}_{dx_}*iHat + "
                      f"{self.source_array}_{dy_}*jHat + "
                      f"{self.source_array}_{dz_}*kHat")
    calc2.UpdatePipeline()

    # reference bounds (un-warped)
    xmin, xmax, ymin, ymax, zmin, zmax = calc2.GetDataInformation().GetBounds()
    n1, n2 = self.resolution
    if self.plane == "xy":
      z = 0.5*(zmin+zmax)
      dims = [n1, n2, 1]
      bounds = [xmin, xmax, ymin, ymax, z, z]
      ds = max((xmax-xmin)/max(n1-1,1), (ymax-ymin)/max(n2-1,1))
    elif self.plane == "xz":
      y = 0.5*(ymin+ymax)
      dims = [n1, 1, n2]
      bounds = [xmin, xmax, y, y, zmin, zmax]
      ds = max((xmax-xmin)/max(n1-1,1), (zmax-zmin)/max(n2-1,1))
    elif self.plane == "yz":
      x = 0.5*(xmin+xmax)
      dims = [1, n1, n2]
      bounds = [x, x, ymin, ymax, zmin, zmax]
      ds = max((ymax-ymin)/max(n1-1,1), (zmax-zmin)/max(n2-1,1))
    else:
      raise ValueError(f"unknown plane: {self.plane}")

    grid = pv.ResampleToImage(Input=calc2)
    grid.UseInputBounds = 0
    grid.SamplingDimensions = dims
    grid.SamplingBounds = bounds
    grid.UpdatePipeline()

    # interpolate rotation AND displacement onto the grid
    interp = pv.PointVolumeInterpolator(Input=calc2, Source=grid)
    interp.Kernel = "GaussianKernel"
    interp.Locator = "Static Point Locator"
    radius = self.kernel_radius if self.kernel_radius is not None else 2*ds
    interp.Kernel.Radius = radius
    interp.UpdatePipeline()

    # warp the grid by the interpolated displacement (matches main Warp)
    warped = pv.WarpByVector(Input=interp)
    warped.Vectors = ["POINTS", "displacement"]
    warped.ScaleFactor = self.warp_scale
    warped.UpdatePipeline()

    # fixed vertical offset on top
    if self.offset_z != 0.0:
      off = pv.Calculator(Input=warped)
      off.AttributeType = "Point Data"
      off.CoordinateResults = 1
      off.Function = f"coords + {self.offset_z}*kHat"
      off.UpdatePipeline()
      glyph_input = off
    else:
      glyph_input = warped

    glyph = pv.Glyph(Input=glyph_input, GlyphType="Arrow")
    glyph.OrientationArray = ["POINTS", "rotation"]
    glyph.ScaleArray = ["POINTS", "rotation"]
    glyph.ScaleFactor = ds / (2. * math.pi) * self.scale
    glyph.GlyphMode = "All Points"
    glyph.UpdatePipeline()

    display = pv.Show(glyph, rview)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb
    return pipe

# --- Runner ------------------------------------------------------------------

def netvis(path, ops=(SolidColor("white")), bg="black", view="iso", resolution=(1000, 1000),
           axis=True, output=None, show=True, duration=5, fps=30):
  if not os.path.isfile(path):
    sys.exit(f"error: file not found: {path}")

  reader = pv.VTKHDFReader(FileName=[path])
  times = reader.TimestepValues

  pipe = reader
  pipe.UpdatePipeline()
  pipe = pv.ExtractSurface(Input=pipe)
  pipe.UpdatePipeline()
  rview = pv.GetActiveViewOrCreate("RenderView")

  for op in ops:
    pipe = op.apply(pipe, rview)

  rview.ViewSize = list(resolution)
  rview.Background = list(to_rgb(bg))
  rview.UseColorPaletteForBackground = 0
  rview.OrientationAxesVisibility = int(axis)

  view.orient(rview)
  pv.ResetCamera()
  pv.Render()

  if show:
    scene = pv.GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.PlayMode = "Snap To TimeSteps"
    scene.NumberOfFrames = len(times)
    if len(times) > 0:
      scene.FramesPerTimestep = max(1, round(duration * fps / len(times)))
    else:
      scene.FramesPerTimestep = 1
    scene.Play()
    pv.Interact()
  if output:
    if len(times) > 0:
      print("saving animation...")
      scene.FramesPerTimestep = 1
      pv.SaveAnimation(output, rview, FrameRate=fps)
    else:
      pv.SaveScreenshot(output, rview, TransparentBackground=1)


def netvis_overlay(path, ops_frames, times, bg="black", view=None,
                   resolution=(1000, 1000), axis=True, output=None, show=True,
                   reference=None):
  if not os.path.isfile(path):
    sys.exit(f"error: file not found: {path}")

  reader = pv.VTKHDFReader(FileName=[path])
  available = list(reader.TimestepValues)
  rview = pv.GetActiveViewOrCreate("RenderView")

  for t, ops in zip(times, ops_frames):
    t_snap = min(available, key=lambda x: abs(x - t)) if available else t
    idx = available.index(t_snap)

    extract = pv.ExtractTimeSteps(Input=reader)
    extract.TimeStepIndices = [idx]
    extract.UpdatePipeline()

    # shift this branch's single timestep to t=0
    shift = pv.TemporalShiftScale(Input=extract)
    shift.PreShift = -t_snap
    shift.UpdatePipeline()

    pipe = pv.ExtractSurface(Input=shift)
    pipe.UpdatePipeline()
    for op in ops:
        pipe = op.apply(pipe, rview)

  rview.ViewTime = 0.0
  rview.ViewSize = list(resolution)
  rview.Background = list(to_rgb(bg))
  rview.UseColorPaletteForBackground = 0
  rview.OrientationAxesVisibility = int(axis)

  view.orient(rview)
  pv.ResetCamera()
  pv.Render()

  if show:
    pv.Interact()
  if output:
    pv.SaveScreenshot(output, rview, TransparentBackground=1)

# --- CLI ---------------------------------------------------------------------

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--fg", default="white")
  p.add_argument("--bg", default="black")
  p.add_argument("--color-by", default=None, help="color by array 'name' or 'name:N'")
  p.add_argument("--color-invert", action="store_true")
  p.add_argument("--color-categories", help="color categories e.g. '1-4,7'", default="")
  p.add_argument("--warp-by", default="values:6,7,8",
                 help="warp by")
  p.add_argument("--warp-scale", type=float, default=1.,
                 help="warp scale")
  p.add_argument("--arrows", type=float, default=1.,
                 help="place arrows (components 6-8 of 'values')")
  p.add_argument("--arrows-offset", type=float, default=0.,
                 help="arrows offset z")
  p.add_argument("--view", choices=list(View.VIEWS), default="top")
  p.add_argument("--resolution", default="1000x1000")
  p.add_argument("-o", "--output", default=None, help="optional screenshot path")
  p.add_argument("--show", type=int, default=1, help="skip Interact()")
  p.add_argument("--axis", type=int, default=1, help="display orientation axis")
  p.add_argument("-r", "--tubes-radius", type=float, default=20,
                 help="tube radius")
  p.add_argument("--tubes-sides", type=int, default=4)
  p.add_argument("--ref", type=int, default=1, help="show reference outline")
  p.add_argument("--ref-opacity", type=float, default=1.)
  p.add_argument("--fps", type=int, default=30, help="target fps")
  p.add_argument("--duration", type=float, default=5, help="target duration of animation")
  p.add_argument("--frames", default=None,
               help="overlay mode: comma-separated times, e.g. '0.0,0.3,0.6'")
  p.add_argument("--frame-colors", default=None,
               help="comma-separated colors, one per frame (default: cycle fg)")

  args = p.parse_args()

  assert("x" in args.resolution)
  resolution = tuple(map(int, args.resolution.split("x")))
  assert(len(resolution) == 2)

  def ops_factory(fg):
    ops = []
    # ref must go before warp
    if args.ref:
      ops.append(Reference(color=args.fg, opacity=args.ref_opacity))
    if args.arrows != 0.:
      ops.append(CoarseArrows(scale=args.arrows, warp_scale=args.warp_scale, offset_z=args.arrows_offset))
    if args.warp_by.lower() != "none":
      temp = args.warp_by.split(":")
      array = temp[0]
      comps = [int(x) for x in temp[1].split(",")]
      if len(temp) > 2:
        normal = temp[2]
      else:
        normal = None
      ops.append(Warp(scale=args.warp_scale, components=comps, source_array=array, normal=normal))
    if args.tubes_radius != 0.:
      ops.append(Tubes(radius=args.tubes_radius, sides=args.tubes_sides))
    if args.color_by:
      ops.append(ArrayColor(args.color_by, fg=fg, invert=args.color_invert, categories=args.color_categories))
    else:
      ops.append(SolidColor(fg))
    return ops

  if args.frames:
    times = [float(s) for s in args.frames.split(",")]
    if args.frame_colors:
      try:
        cmap = plt.get_cmap(args.frame_colors)
        colors = [cmap(i / max(len(times) - 1, 1)) for i in range(len(times))]
      except ValueError:
        colors = args.frame_colors.split(",")
    else:
      colors = [args.fg] * len(times)

    ops = [ops_factory(fg) for fg in colors]

    netvis_overlay(args.input, ops, times, bg=args.bg,
                   view=View(args.view), axis=args.axis, resolution=resolution,
                   output=args.output, show=args.show)
  else:
    ops = ops_factory(args.fg)

    netvis(args.input, ops=ops, bg=args.bg, view=View(args.view), axis=args.axis,
         resolution=resolution, output=args.output, show=args.show, duration=args.duration, fps=args.fps)

  cam = pv.GetActiveCamera()
  pos   = np.array(cam.GetPosition())
  focal = np.array(cam.GetFocalPoint())
  up    = np.array(cam.GetViewUp())

  rot   = pos - focal
  rot  /= np.linalg.norm(rot)
  up   /= np.linalg.norm(up)

  print("view:")
  print(f'"custom": {{"position": {tuple(rot.round(4).tolist())}, ' f'"focal": (0, 0, 0), "up": {tuple(up.round(4).tolist())}}},')
