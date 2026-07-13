#!/usr/bin/env pvpython

"""
netvis.py — render and animate VTKHDF beam/fiber networks via ParaView.

Builds a ParaView pipeline from a configurable list of ops (Warp, Tubes,
CoarseArrows, Q1Mesh, Reference, SolidColor, ArrayColor) and either opens an
interactive view, dumps a screenshot, or saves a transient animation across
the file's time steps. Single-file mode renders one pipeline over all time
steps; overlay mode (--frames) renders several time steps simultaneously in
the same view, each as its own pipeline branch shifted to t=0.

Typical use is to visualize HDG basis functions / eigenmodes / displacement
columns produced by netcoarse.py, with the 6/7/8 components of the 'values'
array interpreted as displacement and 9/10/11 as rotation (Timoshenko DoFs).

Examples
--------
    netvis.py field.vtkhdf -r 20 --view iso --warp-scale 1e3
    netvis.py field.vtkhdf --color-by values:0 --warp-by none --arrows 0
    netvis.py field.vtkhdf --frames 0.0,0.3,0.6 --frame-colors viridis

Inputs
------
VTKHDF UnstructuredGrid of VTK_LINE cells, typically with a multi-component
'values' PointData array and one or more time steps.

Output
------
Interactive RenderView, single screenshot, or animation file, depending on
--show / --output and whether the file has time steps.

Author
------
Joseph Holten, KIT, 2026.
"""

import argparse
import sys
import os
# use every core for VTK's SMP-parallel filters (must be set before VTK loads);
# respects an externally set value
os.environ.setdefault("VTK_SMP_MAX_THREADS", str(os.cpu_count() or 1))
import paraview.simple as pv
try:
  from matplotlib.colors import to_rgb
except ImportError:  # ParaView-bundled pythons ship without matplotlib
  _NAMED = {"white": (1., 1., 1.), "black": (0., 0., 0.), "red": (1., 0., 0.),
            "green": (0., .5, 0.), "blue": (0., 0., 1.), "gray": (.5, .5, .5),
            "grey": (.5, .5, .5)}
  def to_rgb(c):
    if isinstance(c, (tuple, list)):
      return tuple(float(x) for x in c[:3])
    if isinstance(c, str) and c.startswith("#") and len(c) == 7:
      return tuple(int(c[i:i + 2], 16) / 255. for i in (1, 3, 5))
    if isinstance(c, str) and c.lower() in _NAMED:
      return _NAMED[c.lower()]
    raise ValueError(f"color '{c}' needs matplotlib (or use #rrggbb)")
from colorsys import hsv_to_rgb
import math
import numpy as np
import builtins

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
        for i in range(info.GetNumberOfArrays()):
            arr = info.GetArray(i)
            if arr.GetName() == name:
                return assoc, arr
    return None, None

# --- Pipeline ops ------------------------------------------------------------
# Each op is a callable: pipe -> pipe. May skip and return input unchanged.

class Warp:
  def __init__(self, components=(6, 7, 8), scale=1.0, source_array="values", normal=None):
    self.components = components
    self.scale = scale
    self.source_array = source_array
    self.normal = normal

  def apply(self, pipe, view):
    assoc, arr = find_array(pipe, self.source_array)
    if assoc != "POINTS":
      print(f"warning: Warp: '{self.source_array}' not found, skipping",
            file=sys.stderr)
      return pipe
    ncomp = arr.GetNumberOfComponents()
    needed = builtins.max(self.components) + 1
    if ncomp < needed:
        print(f"warning: Warp: '{self.source_array}' has {ncomp} "
              f"components, need {needed}, skipping", file=sys.stderr)
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

class Surface:
  """Triangulate the (already warped) points into a filled 2D surface.

  The right rendering for planar networks whose edges are denser than pixels (fine
  grids alias into moire as lines/tubes): only the surface texture and the global
  displacement trend are visible anyway, so fill them in. Merges the per-edge
  duplicate endpoints first (the plot writes each node once per adjacent edge), then
  Delaunay-triangulates the xy projection. alpha > 0 bounds the triangle
  circumradius so concave regions/holes are not bridged (use ~2x the mean edge
  length for disordered nets; 0 = fill to the convex hull, fine for grids).
  """
  def __init__(self, alpha=0.0):
    self.alpha = alpha

  def apply(self, pipe, view):
    clean = pv.Clean(Input=pipe)
    clean.UpdatePipeline()
    tri = pv.Delaunay2D(Input=clean)
    if self.alpha > 0:
      tri.Alpha = self.alpha
    tri.UpdatePipeline()
    return tri

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

class Beams:
  """Build hollow rectangular beams (4 side quads per edge) from CellData normals/widths.

  Reads two 3-vector normals and two scalar widths per cell from `source_array`
  (CellData). For each VTK_LINE cell with endpoints p0, p1, places 4 corners
  around each endpoint at p_i +/- (w1/2) n1 +/- (w2/2) n2, then emits the 4
  side quads connecting them. Each quad inherits the source edge's CellData
  (4 copies per edge) so per-beam coloring via --color-by works directly.

  Run AFTER Warp: the centerline endpoints are already displaced by the
  per-node displacement, so all 4 corners of a section share that displacement
  by construction. Fixed cross-section orientation (does not rotate with the
  beam's rotation DoF).
  """
  def __init__(self, source_array="properties",
               n1_cols=(7,8,9), n2_cols=(10,11,12),
               w1_col=13, w2_col=14, scale=1.0,
               skip_col=15, skip_value=-1.0,
               rot_array="values", rot_cols=(9,10,11), rot_scale=1.0):
    self.source_array = source_array
    self.n1_cols = tuple(n1_cols)
    self.n2_cols = tuple(n2_cols)
    self.w1_col = w1_col
    self.w2_col = w2_col
    self.scale = scale
    self.skip_col = skip_col         # column in source_array used to drop edges; None disables
    self.skip_value = skip_value     # edges where source_array[:, skip_col] == skip_value are dropped
    self.rot_array = rot_array       # PointData array holding rotation vectors; "" disables
    self.rot_cols = tuple(rot_cols)  # 3 columns of rot_array taken as the rotation vector
    self.rot_scale = rot_scale       # multiplier applied to rotation vectors (visualization gain)

  def apply(self, pipe, rview):
    pf = pv.ProgrammableFilter(Input=pipe)
    pf.OutputDataSetType = "vtkPolyData"
    pf.Script = f"""
import numpy as np
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkCellArray
from vtkmodules.util.vtkConstants import VTK_ID_TYPE
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy

vin  = self.GetInputDataObject(0, 0)
vout = self.GetOutputDataObject(0)
vout.Initialize()   # discard verts/lines/polys/strips and arrays shallow-copied from the input

prop = vin.GetCellData().GetArray("{self.source_array}")
if prop is None:
    raise RuntimeError("Beams: CellData array '{self.source_array}' not found")
prop = vtk_to_numpy(prop)

pts_in = vtk_to_numpy(vin.GetPoints().GetData())
nc = vin.GetNumberOfCells()

n1 = prop[:, [{self.n1_cols[0]}, {self.n1_cols[1]}, {self.n1_cols[2]}]]
n2 = prop[:, [{self.n2_cols[0]}, {self.n2_cols[1]}, {self.n2_cols[2]}]]
w1 = prop[:, {self.w1_col}] * {self.scale}
w2 = prop[:, {self.w2_col}] * {self.scale}

# Vectorized read of line endpoints via the polydata Lines connectivity.
# For a polydata of pure VTK_LINE cells, GetConnectivityArray() is a flat
# (2*nc,) int array of [p0,p1, p0,p1, ...].
lines = vin.GetLines()
conn  = vtk_to_numpy(lines.GetConnectivityArray()).reshape(nc, 2)
cell_pids = conn.astype(np.int64)

# Drop cells where source_array[:, skip_col] == skip_value (e.g. fiber_id == -1 connectors).
keep = np.ones(nc, dtype=bool)
skip_col = {self.skip_col!r}
if skip_col is not None:
    keep = prop[:, skip_col] != {self.skip_value}
keep_idx  = np.flatnonzero(keep)              # original cell indices kept
n1        = n1[keep]
n2        = n2[keep]
w1        = w1[keep]
w2        = w2[keep]
cell_pids = cell_pids[keep]
nc        = int(keep.sum())

p0 = pts_in[cell_pids[:, 0]]
p1 = pts_in[cell_pids[:, 1]]

# Per-endpoint rotated normals via Rodrigues' formula on a per-node rotation vector.
# v_rot = v cos t + (k x v) sin t + k (k . v) (1 - cos t),  with k = r/|r|, t = |r|.
def _rodrigues(r, v):
    t = np.linalg.norm(r, axis=1, keepdims=True)
    safe = np.where(t > 0, t, 1.0)
    k = r / safe
    kxv = np.cross(k, v)
    kdv = np.sum(k * v, axis=1, keepdims=True)
    c, s = np.cos(t), np.sin(t)
    return v * c + kxv * s + k * kdv * (1.0 - c)

rot_name = "{self.rot_array}"
rot = None
if rot_name:
    ra = vin.GetPointData().GetArray(rot_name)
    if ra is not None:
        rv = vtk_to_numpy(ra)
        rot = rv[:, [{self.rot_cols[0]}, {self.rot_cols[1]}, {self.rot_cols[2]}]] * {self.rot_scale}
    if ra is None:
        print("Beams: rotate: '{self.rot_array}' not found, skipping")

if rot is not None:
    r0 = rot[cell_pids[:, 0]]
    r1 = rot[cell_pids[:, 1]]
    n1_0 = _rodrigues(r0, n1)
    n2_0 = _rodrigues(r0, n2)
    n1_1 = _rodrigues(r1, n1)
    n2_1 = _rodrigues(r1, n2)
else:
    n1_0 = n1_1 = n1
    n2_0 = n2_1 = n2

h1_0 = 0.5 * w1[:, None] * n1_0
h2_0 = 0.5 * w2[:, None] * n2_0
h1_1 = 0.5 * w1[:, None] * n1_1
h2_1 = 0.5 * w2[:, None] * n2_1

# 4 corners going around: (+,+) (-,+) (-,-) (+,-)
s1 = np.array([+1, -1, -1, +1])[None, :, None]
s2 = np.array([+1, +1, -1, -1])[None, :, None]
c0 = p0[:, None, :] + s1*h1_0[:, None, :] + s2*h2_0[:, None, :]   # (nc,4,3)
c1 = p1[:, None, :] + s1*h1_1[:, None, :] + s2*h2_1[:, None, :]
# layout: per-cell corners [0..3] at p0, [4..7] at p1
corners = np.concatenate([c0, c1], axis=1).reshape(-1, 3)

# 4 side quads per cell, winding consistently around the section
k  = np.arange(4)
kn = (k + 1) % 4
base = (np.arange(nc) * 8)[:, None]
quads = np.stack([
    base + k[None, :],
    base + kn[None, :],
    base + 4 + kn[None, :],
    base + 4 + k[None, :],
], axis=2).reshape(-1, 4).astype(np.int64)
nq = quads.shape[0]

vpts = vtkPoints()
vpts.SetData(numpy_to_vtk(np.ascontiguousarray(corners, dtype=np.float64), deep=1))
vout.SetPoints(vpts)

offsets      = np.arange(0, (nq + 1) * 4, 4, dtype=np.int64)
connectivity = np.ascontiguousarray(quads.ravel(), dtype=np.int64)
ca = vtkCellArray()
ca.SetData(numpy_to_vtk(offsets,      deep=1, array_type=VTK_ID_TYPE),
           numpy_to_vtk(connectivity, deep=1, array_type=VTK_ID_TYPE))
vout.SetPolys(ca)

# CellData: each output quad inherits from its source edge (4 copies per edge)
src_cid = np.repeat(keep_idx, 4)
in_cd, out_cd = vin.GetCellData(), vout.GetCellData()
for ai in range(in_cd.GetNumberOfArrays()):
    a = in_cd.GetArray(ai)
    if a is None: continue
    propagated = vtk_to_numpy(a)[src_cid]
    va = numpy_to_vtk(np.ascontiguousarray(propagated), deep=1)
    va.SetName(a.GetName())
    out_cd.AddArray(va)

# PointData: each corner inherits from its source centerline node
src_pid = np.empty((nc, 8), dtype=np.int64)
src_pid[:, :4] = cell_pids[:, 0:1]
src_pid[:, 4:] = cell_pids[:, 1:2]
src_pid = src_pid.ravel()
in_pd, out_pd = vin.GetPointData(), vout.GetPointData()
for ai in range(in_pd.GetNumberOfArrays()):
    a = in_pd.GetArray(ai)
    if a is None: continue
    propagated = vtk_to_numpy(a)[src_pid]
    va = numpy_to_vtk(np.ascontiguousarray(propagated), deep=1)
    va.SetName(a.GetName())
    out_pd.AddArray(va)
"""
    pf.UpdatePipeline()
    return pf

class Q1Mesh:
  def __init__(self, dims, color="magenta", opacity=1.0, line_width=1.0, scale=(1,1,1), eps=0.01, offset=(0,0,0)):
    self.dims = dims
    self.color = color
    self.opacity = opacity
    self.line_width = line_width
    self.scale = scale
    self.eps = eps
    self.offset = offset

  def apply(self, pipe, rview):
    xmin, xmax, ymin, ymax, zmin, zmax = pipe.GetDataInformation().GetBounds()
    nx, ny, nz = self.dims

    # rescale bounds by self.scale
    # extent in any direction should be at least self.eps * max extent
    cx, cy, cz = 0.5*(xmin+xmax), 0.5*(ymin+ymax), 0.5*(zmin+zmax)
    hx, hy, hz = 0.5*(xmax-xmin), 0.5*(ymax-ymin), 0.5*(zmax-zmin)
    eps = self.eps * builtins.max(hx, hy, hz)
    hx = hx or eps
    hy = hy or eps
    hz = hz or eps
    sx, sy, sz = self.scale
    ox, oy, oz = self.offset
    xmin, xmax = cx - sx*hx + ox, cx + sx*hx + ox
    ymin, ymax = cy - sy*hy + oy, cy + sy*hy + oy
    zmin, zmax = cz - sz*hz + oz, cz + sz*hz + oz

    # Wavelet produces an ImageData with the given extent, centered at Center
    src = pv.Wavelet()
    src.WholeExtent = [0, nx, 0, ny, 0, builtins.max(nz, 0)]
    # Place origin at 0 by setting Center to half-extent
    src.Center = [nx / 2.0, ny / 2.0, builtins.max(nz, 0) / 2.0]
    src.UpdatePipeline()

    # Scale + translate to target bounds
    sx = (xmax - xmin) / builtins.max(nx, 1)
    sy = (ymax - ymin) / builtins.max(ny, 1)
    sz = (zmax - zmin) / builtins.max(nz, 1) if nz > 0 and zmax > zmin else 1.0

    tf = pv.Transform(Input=src)
    tf.Transform = "Transform"
    tf.Transform.Scale = [sx, sy, sz]
    tf.Transform.Translate = [xmin, ymin, zmin]
    tf.UpdatePipeline()

    display = pv.Show(tf, rview)
    display.Representation = "Wireframe"
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb
    display.Opacity = self.opacity
    display.LineWidth = self.line_width
    return pipe

class Reference:
  def __init__(self, color="white", opacity=1.0, line_width=1.0):
    self.color = color
    self.opacity = opacity
    self.line_width = line_width

  def apply(self, pipe, view):
    outline = pv.Outline(Input=pipe)
    outline.UpdatePipeline()
    display = pv.Show(outline, view)
    rgb = list(to_rgb(self.color))
    display.AmbientColor = rgb
    display.DiffuseColor = rgb
    display.Opacity = self.opacity
    display.LineWidth = self.line_width
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


# (display, ctf) pairs whose transfer function must be rescaled per animation frame
# (ArrayColor rescale="frame"); consumed by the frame loop in netvis()
FRAME_RESCALE = []

class ArrayColor:
  def __init__(self, spec, fg="white", invert=False, categories="", bg="black",
               rescale="time"):
    """'name' or 'name:N' -> (name, component_or_None)."""
    self.invert = invert
    self.fg = fg
    self.bg = bg
    self.rescale = rescale

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
    assoc, arr = find_array(pipe, self.name)
    if assoc is None:
      print(f"warning: '{self.name}' not found", file=sys.stderr)
      return
    if self.comp is None:
      target = (assoc, self.name)
    else:
      target = (assoc, self.name, self.comp)
    pv.ColorBy(display, target)
    ctf = pv.GetColorTransferFunction(self.name)
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
      if self.rescale == "frame":
        # rescale to the CURRENT step now and register for the animation loop:
        # as the wave spreads, the amplitude drops and a global range hides the
        # front in later frames; the trade-off is that the colorbar changes
        # meaning between frames.
        display.RescaleTransferFunctionToDataRange(False, True)
        FRAME_RESCALE.append((display, ctf, self.name))
      else:
        # "Rescale to Data Range Over All Timesteps": ParaView sweeps every
        # registered time step in C++ and rescales the LUT to the global range
        # of the colored component. The per-display array info only sees the
        # current step, so without this the colors jump every frame.
        display.RescaleTransferFunctionToDataRangeOverTime()
      # symmetrize the diverging scale around 0 (RGBPoints is a flat
      # [scalar, r, g, b, ...] list, so [0] / [-4] are the rescaled min / max)
      pts = ctf.RGBPoints
      M = builtins.max(abs(pts[0]), abs(pts[-4]))
      ctf.ApplyPreset("Cool to Warm", True)
      ctf.RescaleTransferFunction(-M, M)
      if self.invert:
        ctf.InvertTransferFunction()
    display.SetScalarBarVisibility(pv.GetActiveView(), True)
    # black scalebar text on light backgrounds (and transparent stills viewed on
    # white); the default white text is invisible there
    lum = sum(to_rgb(self.bg)[:3]) / 3.
    txt = [0., 0., 0.] if lum > 0.5 else [1., 1., 1.]
    sb = pv.GetScalarBar(ctf, rview)
    sb.TitleColor = txt
    sb.LabelColor = txt

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
    assoc, arr = find_array(pipe, self.source_array)
    if assoc != "POINTS":
      print(f"warning: CoarseArrows: '{self.source_array}' not found, skipping",
        file=sys.stderr)
      return pipe
    ncomp = arr.GetNumberOfComponents()
    needed = builtins.max(*self.disp_components, *self.rot_components) + 1
    if ncomp < needed:
        print(f"warning: CoarseArrows: '{self.source_array}' has {ncomp} "
              f"components, need {needed}, skipping", file=sys.stderr)
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
      z = 1.2 * zmax
      dims = [n1, n2, 1]
      bounds = [xmin, xmax, ymin, ymax, z, z]
      ds = builtins.max((xmax-xmin)/builtins.max(n1-1,1), (ymax-ymin)/builtins.max(n2-1,1))
    elif self.plane == "xz":
      y = 1.2 * ymax
      dims = [n1, 1, n2]
      bounds = [xmin, xmax, y, y, zmin, zmax]
      ds = builtins.max((xmax-xmin)/builtins.max(n1-1,1), (zmax-zmin)/builtins.max(n2-1,1))
    elif self.plane == "yz":
      x = 1.2 * xmax
      dims = [1, n1, n2]
      bounds = [x, x, ymin, ymax, zmin, zmax]
      ds = builtins.max((ymax-ymin)/builtins.max(n1-1,1), (zmax-zmin)/builtins.max(n2-1,1))
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
  """Render `path` with a single pipeline, optionally animating across time steps.

  Applies `ops` in order to a surface-extracted reader, configures the
  RenderView (background, size, orientation axis, camera via `view`), then:
    - if `show`: plays the animation in an interactive window;
    - if `output`: saves an animation (if time steps exist) or a single
      screenshot. Animation length is approximately `duration` seconds at
      `fps`, snapped to the file's time steps.
  """
  if not os.path.isfile(path):
    sys.exit(f"error: file not found: {path}")

  reader = pv.VTKHDFReader(FileName=[path])
  times = reader.TimestepValues

  pipe = reader
  pipe.UpdatePipeline()
  # under MPI (mpirun -np N pvbatch netvis.py ...) the reader loads a single part;
  # RedistributeDataSet partitions it once so every downstream filter runs N-way
  # data-parallel with IceT compositing the render. Serial pvpython is unaffected.
  # Caveat: --surface triangulates per rank, so partition-boundary seams can appear.
  from paraview import servermanager
  nranks = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
  if nranks > 1:
    print(f"MPI: redistributing data over {nranks} ranks")
    pipe = pv.RedistributeDataSet(Input=pipe)
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
      scene.FramesPerTimestep = builtins.max(1, int(round(duration * fps / len(times))))
    else:
      scene.FramesPerTimestep = 1
    scene.Play()
    pv.Interact()
  if output:
    if len(times) > 0:
      # headless (--show 0): the scene setup above was skipped, do it here
      scene = pv.GetAnimationScene()
      scene.UpdateAnimationUsingDataTimeSteps()
      scene.PlayMode = "Snap To TimeSteps"
      scene.NumberOfFrames = len(times)
      scene.FramesPerTimestep = 1
      # one frame per time step; stretch playback to ~duration via the file fps
      # (30 fps for 21 steps = 0.7 s of video, unwatchable)
      rate = builtins.max(1, int(round(len(times) / duration)))
      if FRAME_RESCALE:
        # dynamic LUT inside SaveAnimation: a PythonAnimationCue ticks once per
        # frame, rescaling each registered transfer function to the CURRENT step's
        # range re-symmetrized about 0. (The built-in AutomaticRescaleRangeMode
        # "Clamp and update every timestep" also rescales per frame but cannot keep
        # the diverging map centered on zero.) The extra Render() in tick forces the
        # pipeline update so the data range belongs to this frame, not the previous.
        names = sorted({name for _, _, name in FRAME_RESCALE})
        cue = pv.PythonAnimationCue()
        cue.Script = f'''
def start_cue(self): pass

def tick(self):
    import paraview.simple as pv
    pv.Render()
    view = pv.GetActiveView()
    for rep in view.Representations:
        an = getattr(rep, "ColorArrayName", None)
        if an is None or an[1] not in {names!r} or not getattr(rep, "Visibility", 0):
            continue
        rep.RescaleTransferFunctionToDataRange(False, True)
        ctf = pv.GetColorTransferFunction(an[1])
        pts = ctf.RGBPoints
        M = max(abs(pts[0]), abs(pts[-4]), 1e-300)   # t=0 is all-zero: avoid a 0-width range
        ctf.RescaleTransferFunction(-M, M)

def end_cue(self): pass
'''
        scene.Cues.append(cue)
        print("saving animation (per-frame color rescale via animation cue)...")
        pv.SaveAnimation(output, rview, FrameRate=rate)
        scene.Cues.remove(cue)
      else:
        print("saving animation...")
        pv.SaveAnimation(output, rview, FrameRate=rate)
    else:
      pv.SaveScreenshot(output, rview, TransparentBackground=1)


def netvis_overlay(path, ops_frames, times, bg="black", view=None,
                   resolution=(1000, 1000), axis=True, output=None, show=True,
                   reference=None):
  """Render several time steps of `path` simultaneously, each with its own ops.

  For each (t, ops) pair, extracts the time step nearest `t`, shifts it to
  t=0 via TemporalShiftScale so all branches coexist in one frame, and
  applies `ops` to that branch. Useful for showing mode shapes side-by-side
  or comparing snapshots of a transient field. Output is a single
  screenshot (no animation in overlay mode).
  """
  if not os.path.isfile(path):
    sys.exit(f"error: file not found: {path}")

  reader = pv.VTKHDFReader(FileName=[path])
  available = list(reader.TimestepValues)
  rview = pv.GetActiveViewOrCreate("RenderView")

  for t, ops in zip(times, ops_frames):
    t_snap = builtins.min(available, key=lambda x: abs(x - t)) if available else t
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
  p = argparse.ArgumentParser(
    description="Render VTKHDF beam/fiber networks via ParaView: tubes, "
                "warping by displacement, coarse rotation arrows, optional "
                "Q1 background mesh, and multi-frame overlay mode.")
  p.add_argument("input", help="path to .vtkhdf file")
  p.add_argument("--fg", default="white", help="foreground color (solid coloring + reference outline)")
  p.add_argument("--bg", default="black", help="background color")
  p.add_argument("--color-by", default=None, help="color by array 'name' or 'name:N'")
  p.add_argument("--color-invert", action="store_true", help="invert the Cool-to-Warm transfer function")
  p.add_argument("--color-rescale", choices=["time", "frame"], default="time",
                 help="color range over all timesteps (default; comparable frames) or "
                      "per frame (keeps the decaying wave front visible in animations, "
                      "but the colorbar changes meaning between frames)")
  p.add_argument("--color-categories", default="",
                 help="treat values as categorical, e.g. '1-4,7'; uses HSV-spaced colors")
  p.add_argument("--warp-by", default="values:6,7,8",
                 help="warp spec 'array:c1,c2,c3' (vector) or 'array:c:axis' (scalar*axis); 'none' disables")
  p.add_argument("--warp-scale", type=float, default=1.,
                 help="scale factor applied to the warp vector")
  p.add_argument("--arrows", type=float, default=1.,
                 help="scale for coarse rotation arrows (components 9-11 of 'values'); 0 disables")
  p.add_argument("--arrows-offset", type=float, default=0.,
                 help="fixed vertical offset for arrows after warping")
  p.add_argument("--view", choices=list(View.VIEWS), default="top",
                 help="camera preset; non-'iso' use parallel projection")
  p.add_argument("--resolution", default="1000x1000", help="render size WxH")
  p.add_argument("-o", "--output", default=None, help="optional screenshot path")
  p.add_argument("--show", type=int, default=1, help="if 1, open interactive window (Interact)")
  p.add_argument("--axis", type=int, default=1, help="display orientation axis")
  p.add_argument("-r", "--tubes-radius", type=float, default=20,
                 help="tube radius for beam rendering; 0 disables tubes")
  p.add_argument("--tubes-sides", type=int, default=4, help="number of polygonal sides per tube")
  p.add_argument("--surface", action="store_true",
                 help="render the network as a filled Delaunay surface of its (warped) "
                      "points instead of tubes/beams; for planar nets finer than pixels "
                      "(fine grids moire as lines)")
  p.add_argument("--surface-alpha", type=float, default=0.,
                 help="Delaunay alpha radius for --surface: triangles above this "
                      "circumradius are dropped (~2x mean edge length keeps holes open "
                      "in disordered nets); 0 fills to the convex hull")
  p.add_argument("--beams", type=int, default=0,
                 help="if 1, render edges as hollow rectangular beams using CellData normals/widths; overrides --tubes-radius")
  p.add_argument("--beams-array", default="properties",
                 help="CellData array holding beam normals and widths")
  p.add_argument("--beams-cols", default="7,8,9:10,11,12:13:14",
                 help="column spec 'n1x,n1y,n1z:n2x,n2y,n2z:w1:w2' (0-indexed)")
  p.add_argument("--beams-scale", type=float, default=1.0,
                 help="scale factor applied to beam widths")
  p.add_argument("--beams-skip", default="15=-1",
                 help="drop edges where source_array[:, COL] == VAL, as 'COL=VAL'; '' disables")
  p.add_argument("--beams-rotate", default="values:9,10,11",
                 help="rotate beam cross-sections via Rodrigues using PointData 'ARRAY:c1,c2,c3'; '' disables")
  p.add_argument("--beams-rot-scale", type=float, default=1.0,
                 help="scale rotation vectors before applying")
  p.add_argument("--ref", type=int, default=1, help="show reference outline")
  p.add_argument("--ref-opacity", type=float, default=1., help="opacity of reference outline")
  p.add_argument("--fps", type=int, default=30, help="target fps")
  p.add_argument("--duration", type=float, default=5, help="target duration of animation")
  p.add_argument("--frames", default=None,
               help="overlay mode: comma-separated times, e.g. '0.0,0.3,0.6'")
  p.add_argument("--frame-colors", default=None,
               help="matplotlib colormap name or comma-separated colors, one per frame")
  p.add_argument("--q1", default=None,
               help="overlay Q1 cartesian mesh: 'N' (N,N,1), 'NX,NY' (NX,NY,1), or 'NX,NY,NZ'")
  p.add_argument("--q1-color", default="black",
               help="overlay Q1 cartesian mesh color")
  p.add_argument("--q1-scale", default="1,1,1",
               help="overlay Q1 cartesian mesh scale")
  p.add_argument("--q1-offset", default="0,0,0",
               help="overlay Q1 cartesian mesh offset")
  p.add_argument("--line-width", type=float, default=1, help="line width for wiremeshes")

  args = p.parse_args()

  assert("x" in args.resolution)
  resolution = tuple(map(int, args.resolution.split("x")))
  assert(len(resolution) == 2)

  def ops_factory(fg):
    ops = []
    # ref must go before warp
    if args.q1:
      parts = [int(x) for x in args.q1.split(",")]
      if len(parts) == 1:
        dims = (parts[0], parts[0], 1)
      elif len(parts) == 2:
        dims = (parts[0], parts[1], 1)
      elif len(parts) == 3:
        dims = tuple(parts)
      else:
        parser.error("error: --q1 takes 1, 2, or 3 comma-separated ints")
      scale = [float(x) for x in args.q1_scale.split(",")]
      if len(scale) != 3:
        parser.error("--q1-scale takes 3 comma-separated floats")
      offset = [float(x) for x in args.q1_offset.split(",")]
      if len(scale) != 3:
        parser.error("--q1-scale takes 3 comma-separated floats")
      ops.append(Q1Mesh(dims, color=args.q1_color, line_width=args.line_width,
                        offset=offset, scale=scale))
    if args.q1 is None and args.ref:
      ops.append(Reference(color=args.fg, opacity=args.ref_opacity, line_width=args.line_width))
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
    if args.surface:
      ops.append(Surface(alpha=args.surface_alpha))
    elif args.beams:
      parts = args.beams_cols.split(":")
      if len(parts) != 4:
        parser.error("--beams-cols takes 'n1x,n1y,n1z:n2x,n2y,n2z:w1:w2'")
      n1_cols = tuple(int(x) for x in parts[0].split(","))
      n2_cols = tuple(int(x) for x in parts[1].split(","))
      if len(n1_cols) != 3 or len(n2_cols) != 3:
        parser.error("--beams-cols: each normal needs 3 comma-separated cols")
      if args.beams_skip:
        if "=" not in args.beams_skip:
          parser.error("--beams-skip takes 'COL=VAL'")
        skip_col, skip_value = args.beams_skip.split("=", 1)
        skip_col = int(skip_col)
        skip_value = float(skip_value)
      else:
        skip_col, skip_value = None, 0.0
      if args.beams_rotate:
        if ":" not in args.beams_rotate:
          parser.error("--beams-rotate takes 'ARRAY:c1,c2,c3'")
        rot_array, rot_cols_s = args.beams_rotate.split(":", 1)
        rot_cols = tuple(int(x) for x in rot_cols_s.split(","))
        if len(rot_cols) != 3:
          parser.error("--beams-rotate: need 3 comma-separated cols")
      else:
        rot_array, rot_cols = "", (0, 0, 0)
      ops.append(Beams(source_array=args.beams_array,
                       n1_cols=n1_cols, n2_cols=n2_cols,
                       w1_col=int(parts[2]), w2_col=int(parts[3]),
                       scale=args.beams_scale,
                       skip_col=skip_col, skip_value=skip_value,
                       rot_array=rot_array, rot_cols=rot_cols, rot_scale=args.beams_rot_scale))
    elif args.tubes_radius != 0.:
      ops.append(Tubes(radius=args.tubes_radius, sides=args.tubes_sides))
    if args.color_by:
      ops.append(ArrayColor(args.color_by, fg=fg, invert=args.color_invert,
                            categories=args.color_categories, bg=args.bg,
                            rescale=args.color_rescale))
    else:
      ops.append(SolidColor(fg))
    return ops

  if args.frames:
    times = [float(s) for s in args.frames.split(",")]
    if args.frame_colors:
      try:
        import matplotlib.pyplot as plt  # optional: only for --frame-colors colormaps
        cmap = plt.get_cmap(args.frame_colors)
        colors = [cmap(i / builtins.max(len(times) - 1, 1)) for i in range(len(times))]
      except (ValueError, ImportError):
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
