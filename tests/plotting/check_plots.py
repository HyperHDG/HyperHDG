"""Read plots of u = t * x through ParaView, check every point value and render the plots.

    pvbatch tests/plotting/check_plots.py output/plot_exact-2-3.vtu output/plot_vtkhdf-1-3.vtkhdf

The files are those of tests/plotting/plot_exact (VTU, plotted at t = 1) and plot_vtkhdf (VTKHDF,
three time steps): every hypernode is Dirichlet with the value t * x, so the last component of the
point array "values" is u = t * x at every point and every time step (t = 1 for files without
time steps). Each file is rendered at its last time step, coloured by u, to output/plots/<name>.png.
Run with pvbatch; the official headless ParaView builds render without a display.
"""

import os
import sys

from paraview import servermanager
from paraview.simple import (ColorBy, Delete, GetActiveViewOrCreate,
                             GetColorTransferFunction, OpenDataFile, ResetCamera, SaveScreenshot,
                             Show)
from vtkmodules.numpy_interface import dataset_adapter as dsa

failed = 0
os.makedirs("output/plots", exist_ok=True)
for path in sys.argv[1:]:
  reader = OpenDataFile(path)
  reader.UpdatePipelineInformation()
  times = list(getattr(reader, "TimestepValues", None) or []) or [None]
  worst, points = 0., 0
  for t in times:
    reader.UpdatePipeline(t) if t is not None else reader.UpdatePipeline()
    data = dsa.WrapDataObject(servermanager.Fetch(reader))
    values = data.PointData["values"]
    u = values[:, -1] if len(values.shape) > 1 else values
    x = data.Points[:, 0]
    worst = max(worst, float(abs(u - (1. if t is None else t) * x).max()))
    points = len(x)
  ok = points > 0 and worst < 1e-5
  name = os.path.splitext(os.path.basename(path))[0]
  print(f"  {path}: {len(times)} step(s), {points} points, max |u - t x| = {worst:.1e}: "
        f"{'ok' if ok else 'FAILED'}")
  failed += not ok

  view = GetActiveViewOrCreate("RenderView")
  view.ViewSize = [800, 600]
  if times[-1] is not None:
    view.ViewTime = times[-1]
  display = Show(reader, view)
  display.LineWidth = 3.  # one-dimensional hyperedges
  ColorBy(display, ("POINTS", "values"))
  lut = GetColorTransferFunction("values")
  if len(values.shape) > 1:  # colour by u, the last component
    lut.VectorMode = "Component"
    lut.VectorComponent = values.shape[1] - 1
  display.RescaleTransferFunctionToDataRange(False, True)
  display.SetScalarBarVisibility(view, True)
  # The camera is shared by all files, so it is set from scratch: oblique for data in 3D, from
  # above otherwise. ResetCamera keeps the direction and fits the data.
  view.CameraFocalPoint = [0.5, 0.5, 0.5]
  if data.Points[:, 2].max() > data.Points[:, 2].min():
    view.CameraPosition, view.CameraViewUp = [2.4, -1.4, 1.9], [0., 0., 1.]
  else:
    view.CameraPosition, view.CameraViewUp = [0.5, 0.5, 3.], [0., 1., 0.]
  ResetCamera(view)
  SaveScreenshot(f"output/plots/{name}.png", view)
  # The colour map "values" is shared by all files: release it before the next file, whose values
  # may have fewer components (ParaView crashes updating a scalar bar for a missing component).
  display.SetScalarBarVisibility(view, False)
  lut.VectorMode = "Magnitude"
  Delete(reader)

sys.exit(1 if failed else 0)
