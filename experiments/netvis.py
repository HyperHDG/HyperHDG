#!/usr/bin/env pvpython
import argparse
import sys
import os
import paraview.simple as pv
from matplotlib.colors import to_rgb
from colorsys import hsv_to_rgb

class View:
  VIEWS = {
    "top":    {"position": (0, 0, 1), "focal": (0, 0, 0), "up": (0, 1, 0)},
    "bottom": {"position": (0, 0,-1), "focal": (0, 0, 0), "up": (1, 0, 0)},
    "side":   {"position": (0,-1, 0), "focal": (0, 0, 0), "up": (0, 0, 1)},
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
  def __init__(self, components=(6, 7, 8), scale=1.0, source_array="values"):
    self.components = components
    self.scale = scale
    self.source_array = source_array

  def apply(self, pipe, view):
    if find_array(pipe, self.source_array) != "POINTS":
      print(f"warning: Warp: '{self.source_array}' not found, skipping",
            file=sys.stderr)
      return pipe
    cx, cy, cz = self.components
    calc = pv.Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "displacement"
    calc.Function = (f"{self.source_array}_{cx}*iHat + "
                     f"{self.source_array}_{cy}*jHat + "
                     f"{self.source_array}_{cz}*kHat")
    calc.UpdatePipeline()
    pipe = pv.WarpByVector(Input=calc)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = self.scale
    pipe.UpdatePipeline()
    return pipe

class Glyphs:
  def __init__(self, components=(9, 10, 11), scale=1.0, source_array="values",
               mode="All Points", stride=1, color="yellow", offset_z=0):
    self.components = components
    self.scale = scale
    self.source_array = source_array
    self.mode = mode      # "All Points" or "Every Nth Point"
    self.stride = stride
    self.color = color
    self.offset_z = offset_z

  def apply(self, pipe, rview):
    if find_array(pipe, self.source_array) != "POINTS":
      print(f"warning: Glyphs: '{self.source_array}' not found, skipping",
            file=sys.stderr)
      return pipe
    cx, cy, cz = self.components
    calc = pv.Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "rotation"
    calc.Function = (f"{self.source_array}_{cx}*iHat + "
                     f"{self.source_array}_{cy}*jHat + "
                     f"{self.source_array}_{cz}*kHat")
    calc.UpdatePipeline()
    off = pv.Calculator(Input=calc)
    off.AttributeType = "Point Data"
    off.CoordinateResults = 1
    off.Function = f"coords + {self.offset_z}*kHat"
    off.UpdatePipeline()
    glyph = pv.Glyph(Input=off, GlyphType="Arrow")
    glyph.OrientationArray = ["POINTS", "rotation"]
    glyph.ScaleArray = ["POINTS", "rotation"]
    glyph.ScaleFactor = self.scale
    glyph.GlyphMode = self.mode
    if self.mode == "Every Nth Point":
        glyph.Stride = self.stride
    glyph.UpdatePipeline()
    display = pv.Show(glyph, rview)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb
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


# --- CLI ---------------------------------------------------------------------

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--fg", default="white")
  p.add_argument("--bg", default="black")
  p.add_argument("--color-by", default=None, help="color by array 'name' or 'name:N'")
  p.add_argument("--color-invert", action="store_true")
  p.add_argument("--color-categories", help="color categories e.g. '1-4,7'", default="")
  p.add_argument("--warp", type=float, default=1.,
                 help="warp by displacement (components 6-8 of 'values')")
  p.add_argument("--glyphs", type=float, default=1.,
                 help="place glyphs (components 6-8 of 'values')")
  p.add_argument("--glyphs-offset", type=float, default=0.,
                 help="glyph offset z")
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

  args = p.parse_args()

  assert("x" in args.resolution)
  resolution = tuple(map(int, args.resolution.split("x")))
  assert(len(resolution) == 2)

  ops = []
  # ref must go before warp
  if args.ref:
    ops.append(Reference(color=args.fg, opacity=args.ref_opacity))
  if args.warp != 0.:
    ops.append(Warp(scale=args.warp))
  if args.glyphs != 0.:
    ops.append(Glyphs(scale=args.glyphs, offset_z=args.glyphs_offset))
  if args.tubes_radius != 0.:
    ops.append(Tubes(radius=args.tubes_radius, sides=args.tubes_sides))
  if args.color_by:
    ops.append(ArrayColor(args.color_by, fg=args.fg, invert=args.color_invert, categories=args.color_categories))
  else:
    ops.append(SolidColor(args.fg))

  netvis(args.input, ops=ops, bg=args.bg, view=View(args.view), axis=args.axis,
         resolution=resolution, output=args.output, show=args.show, duration=args.duration, fps=args.fps)
