#!/usr/bin/env pvpython

"""
nethist.py — coarse 2D histogram of a VTKHDF beam-network field.

A lightweight cousin of netvis: instead of rendering warped tubes/beams colored
by a 'values' component, it bins edge midpoints onto a regular xy grid, sums the
selected component over the edges falling in each bin, and shows sum/area as a
2D density map (matplotlib). Handy for a quick top-down read of a ripple / mode
shape without the full ParaView render pipeline.

Per VTK_LINE cell the edge value is the mean of its two endpoint values, and the
edge is dropped into the bin containing its midpoint (reference / un-warped
coords). Transient files animate across all time steps with a fixed, symmetric
color scale; single-time files (or --time) render a single frame.

Examples
--------
    nethist.py wave.vtkhdf --color-by values:8 --bins 60 -o hist.mp4
    nethist.py wave.vtkhdf --color-by values:8 --time 0.3 -o hist.png

Author
------
Joseph Holten, KIT, 2026.
"""

import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import paraview.simple as pv
from paraview import servermanager as sm
from vtkmodules.util.numpy_support import vtk_to_numpy


def parse_spec(spec):
  """'name' or 'name:N' -> (name, component)."""
  if ":" in spec:
    name, comp = spec.rsplit(":", 1)
    return name, int(comp)
  return spec, 0


def read_frame(reader, name, comp, t=None):
  """Fetch (edge midpoints, edge values) for one time step.

  Edge value = mean of the selected component over the two endpoints; midpoint =
  mean of the two endpoint positions (reference coords, no warp).
  """
  if t is None:
    reader.UpdatePipeline()
  else:
    reader.UpdatePipeline(t)
  d = sm.Fetch(reader)

  pts = vtk_to_numpy(d.GetPoints().GetData())
  conn = vtk_to_numpy(d.GetCells().GetConnectivityArray()).reshape(-1, 2)
  arr = d.GetPointData().GetArray(name)
  if arr is None:
    sys.exit(f"error: PointData array '{name}' not found")
  val = vtk_to_numpy(arr)
  if val.ndim > 1:
    val = val[:, comp]

  mid = 0.5 * (pts[conn[:, 0]] + pts[conn[:, 1]])
  w = 0.5 * (val[conn[:, 0]] + val[conn[:, 1]])
  return mid[:, :2], w


def histogram(mid, w, bins, bounds):
  """Sum of edge values per bin, divided by bin area (a 2D density)."""
  (xmin, xmax, ymin, ymax) = bounds
  H, xe, ye = np.histogram2d(mid[:, 0], mid[:, 1], bins=bins,
                             range=[[xmin, xmax], [ymin, ymax]], weights=w)
  area = (xe[1] - xe[0]) * (ye[1] - ye[0])
  # transpose so rows are indexed by y for imshow(origin="lower")
  return H.T / area


if __name__ == "__main__":
  p = argparse.ArgumentParser(
    description="Coarse 2D histogram of a VTKHDF beam-network field: bin edge "
                "midpoints, sum a 'values' component per bin, show sum/area.")
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--color-by", default="values:8", help="array 'name' or 'name:N'")
  p.add_argument("--bins", default="60", help="bin count 'N' or 'NX,NY'")
  p.add_argument("--cmap", default="coolwarm", help="matplotlib colormap")
  p.add_argument("--vmax", type=float, default=None,
                 help="symmetric color limit; default = max|density| over all steps")
  p.add_argument("--time", type=float, default=None,
                 help="render a single step nearest this time instead of animating")
  p.add_argument("-o", "--output", default=None, help="output image or video path")
  p.add_argument("--show", type=int, default=1, help="if 1, open an interactive window")
  p.add_argument("--fps", type=int, default=30, help="animation frame rate")
  p.add_argument("--dpi", type=int, default=150, help="output resolution")
  args = p.parse_args()

  if not os.path.isfile(args.input):
    sys.exit(f"error: file not found: {args.input}")

  parts = [int(x) for x in args.bins.split(",")]
  bins = parts if len(parts) == 2 else parts[0]

  name, comp = parse_spec(args.color_by)
  reader = pv.VTKHDFReader(FileName=[args.input])
  reader.UpdatePipeline()
  times = list(reader.TimestepValues) if reader.TimestepValues is not None else []

  # bounds: reference coords are constant in time, so take them from the first frame
  mid0, _ = read_frame(reader, name, comp, times[0] if times else None)
  bounds = (mid0[:, 0].min(), mid0[:, 0].max(), mid0[:, 1].min(), mid0[:, 1].max())
  extent = [bounds[0], bounds[1], bounds[2], bounds[3]]

  # pick the steps to render
  if not times:
    steps = [None]
  elif args.time is not None:
    steps = [min(times, key=lambda x: abs(x - args.time))]
  else:
    steps = times

  # precompute every frame so the color scale can be fixed over time
  frames = []
  for t in steps:
    mid, w = read_frame(reader, name, comp, t)
    frames.append(histogram(mid, w, bins, bounds))

  vmax = args.vmax if args.vmax is not None else max(np.abs(f).max() for f in frames)
  vmin = -vmax

  fig, ax = plt.subplots()
  im = ax.imshow(frames[0], origin="lower", extent=extent, cmap=args.cmap,
                 vmin=vmin, vmax=vmax, aspect="equal", interpolation="nearest")
  fig.colorbar(im, ax=ax, label=f"{args.color_by}")
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  title = ax.set_title(f"{name}:{comp}" + (f"   t = {steps[0]:.3e}" if steps[0] is not None else ""))

  if len(frames) > 1:
    def update(i):
      im.set_data(frames[i])
      title.set_text(f"{name}:{comp}   t = {steps[i]:.3e}")
      return im, title
    anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                   interval=1000 / args.fps, blit=False)
    if args.output:
      print(f"saving animation to {args.output} ...")
      anim.save(args.output, fps=args.fps, dpi=args.dpi)
  elif args.output:
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")

  if args.show:
    plt.show()
