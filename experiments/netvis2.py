#!/usr/bin/env pvpython
import argparse
from paraview.simple import *
from matplotlib.colors import to_rgb

VIEWS = {
  "top":    {"position": (0, 0, 1), "focal": (0, 0, 0), "up": (0, 1, 0)},
  "bottom": {"position": (0, 0,-1), "focal": (0, 0, 0), "up": (1, 0, 0)},
  "side":   {"position": (0,-1, 0), "focal": (0, 0, 0), "up": (0, 0, 1)},
  "iso":    {"position": (1, 1, 1), "focal": (0, 0, 0), "up": (0, 0, 1)},
}

def apply_view(cam, v):
  cam.SetPosition(*v["position"])
  cam.SetFocalPoint(*v["focal"])
  cam.SetViewUp(*v["up"])

def netvis(path, fg="white", bg="black", view="iso", resolution=(1000, 1000), axis=True,
           warp=True, warp_scale=1.0, output=None, show=True):
  pipe = VTKHDFReader(FileName=[path])
  pipe.UpdatePipeline()

  if warp:
    calc = Calculator(Input=pipe)
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "displacement"
    calc.Function = "values_6*iHat + values_7*jHat + values_8*kHat"
    calc.UpdatePipeline()

    pipe = WarpByVector(Input=calc)
    pipe.Vectors = ["POINTS", "displacement"]
    pipe.ScaleFactor = warp_scale
    pipe.UpdatePipeline()

  pipe = Show(pipe, GetActiveViewOrCreate("RenderView"))
  pipe.Representation = "Wireframe"
  pipe.SetScalarColoring(None, 0)
  pipe.AmbientColor = list(to_rgb(fg))
  pipe.DiffuseColor = list(to_rgb(fg))

  pipe = GetActiveView()
  pipe.ViewSize = list(resolution)
  pipe.Background = list(to_rgb(bg))
  pipe.UseColorPaletteForBackground = 0
  pipe.OrientationAxesVisibility = int(axis)
  apply_view(GetActiveCamera(), VIEWS[view])

  ResetCamera()
  Render()
  if show: Interact()
  if output: SaveScreenshot(output, rview, TransparentBackground=1)

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--fg", default="white")
  p.add_argument("--bg", default="black")
  p.add_argument("--view", choices=list(VIEWS), default="iso")
  p.add_argument("--resolution", default="1000x1000")
  p.add_argument("-o", "--output", default=None, help="optional screenshot path")
  p.add_argument("--no-show", action="store_true", help="skip Interact()")
  p.add_argument("--no-axis", action="store_true", help="no axis")
  args = p.parse_args()

  resolution = tuple(map(int, args.resolution.split("x")))
  netvis(args.input, fg=args.fg, bg=args.bg, view=args.view, axis=not args.no_axis,
         resolution=resolution, output=args.output, show=not args.no_show)
