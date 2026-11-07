"""Reference eigenpairs of the thin T-shaped domain by separation of variables and HDG in 2D.

The thin domain is Omega_eps = D_eps x (-1/2, 1/2), where the cross section D_eps (in the plane
z = const) consists of
  plate 1: (2 eps, 1) x (-3/4 eps, 3/4 eps),     plate 2: (-1, -2 eps) x (-3/4 eps, 3/4 eps),
  plate 3: (-3/2 eps, 3/2 eps) x (2 eps, 1),
  junction eps Y_N = (-2 eps, 2 eps) x (-3/4 eps, 3/4 eps) U (-3/2 eps, 3/2 eps) x (-3/4 eps, 2 eps).
The coefficients are kappa = 1, rho = 1 on the plates, and rho = 4 mu / (39 eps) on the junction.
The plate ends x = +-1 and y = 1 are clamped (homogeneous Dirichlet), and all other boundaries
(including both caps z = +-1/2) are natural. Since neither D_eps nor the coefficients depend on z,
the eigenfunctions are u_eps(x, y, z) = w(x, y) cos(k pi (z + 1/2)), k = 0, 1, 2, ..., where
  -div(grad w) + k^2 pi^2 w = lambda rho w   in D_eps.
The spectrum of the thin domain is the union of the spectra of these two-dimensional problems,
which are solved by the HDG method (DiffusionEigsWeighted with reaction c = k^2 pi^2) on tensor
product meshes of D_eps that are graded towards the junction and towards the reentrant corners.
"""

import math
import os

import numpy as np

import common
import limit

PLATE_HALF_THICKNESS = (0.75, 0.75, 1.5)  # d_i / 2 (in units of eps)


def graded(start, end, first, growth=1.5, largest=1. / 16.):
  """Points from start to end with steps first, first * growth, ... (at most largest)."""
  points, step = [start], first
  while points[-1] + step < end - 1e-14:
    points.append(points[-1] + step)
    step = min(step * growth, largest)
  if len(points) > 1 and end - points[-1] < 0.3 * step:
    points[-1] = end
  else:
    points.append(end)
  return points


def corner_grading(corner, left, right, levels, ratio=0.25):
  """Points that are geometrically graded towards corner from both sides (within (left, right))."""
  return [corner - (corner - left) * ratio ** l for l in range(1, levels + 1)] + \
         [corner + (right - corner) * ratio ** l for l in range(1, levels + 1)]


def breakpoints(eps, corner_levels=4, first=None, largest=1. / 16.):
  """Breakpoints of the tensor product mesh in x and y."""
  first = first if first is not None else eps / 2.
  # plates in r-direction: start at 2 eps, contain 4 eps (end of the reflected collar)
  radial = sorted(set([2. * eps + l * eps / 2. for l in range(5)] +
                      graded(4. * eps, 1., first=first, largest=largest)))
  junction_x = [-1.5, -0.75, 0., 0.75, 1.5]
  xs = set(-r for r in radial) | set(radial) | set(eps * x for x in junction_x)
  xs |= set(corner_grading(1.5 * eps, 0.75 * eps, 2. * eps, corner_levels))
  xs |= set(corner_grading(-1.5 * eps, -2. * eps, -0.75 * eps, corner_levels))
  ys = set(eps * y for y in [-0.75, -0.375, 0., 0.375, 0.75, 1.375]) | set(radial)
  ys |= set(corner_grading(0.75 * eps, 0., 2. * eps, corner_levels))
  return np.array(sorted(xs)), np.array(sorted(ys))


def region(x, y, eps):
  """0, 1, 2 for the plates, 3 for the junction, -1 outside of the cross section."""
  if abs(y) < 0.75 * eps:
    if x > 2. * eps:
      return 0
    if x < -2. * eps:
      return 1
    return 3
  if 0. < y and abs(x) < 1.5 * eps:
    return 2 if y > 2. * eps else 3
  return -1


def write_domain(filename, eps, k, corner_levels=4, first=None, largest=1. / 16.):
  xs, ys = breakpoints(eps, corner_levels, first, largest)
  rho_junction = 4. * limit.MU / (39. * eps)
  cells = {}
  for i in range(len(xs) - 1):
    for j in range(len(ys) - 1):
      reg = region(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]), eps)
      if reg >= 0:
        cells[(i, j)] = reg
  points, nodes, edges = {}, {}, []

  def point(i, j):
    return points.setdefault((i, j), len(points))

  def node(key):
    return nodes.setdefault(key, len(nodes))

  def face_type(key):
    kind, i, j = key
    if kind == "v":  # vertical face x = xs[i], ys[j] < y < ys[j + 1]
      neighbors = [(i - 1, j) in cells, (i, j) in cells]
      if all(neighbors):
        return common.INTERIOR
      return common.DIRICHLET if abs(abs(xs[i]) - 1.) < 1e-14 else common.NATURAL
    neighbors = [(i, j - 1) in cells, (i, j) in cells]  # horizontal face y = ys[j]
    if all(neighbors):
      return common.INTERIOR
    return common.DIRICHLET if abs(ys[j] - 1.) < 1e-14 else common.NATURAL

  for (i, j), reg in sorted(cells.items()):
    faces = [("v", i, j), ("v", i + 1, j), ("h", i, j), ("h", i, j + 1)]
    edges.append({
      "points": [point(i, j), point(i + 1, j), point(i, j + 1), point(i + 1, j + 1)],
      "nodes": [node(key) for key in faces],
      "types": [face_type(key) for key in faces],
      "props": [1., 1., rho_junction if reg == 3 else 1., (k * math.pi) ** 2, 0.],
    })
  coords = np.zeros((len(points), 2))
  for (i, j), index in points.items():
    coords[index] = [xs[i], ys[j]]
  common.write_geo(filename, 2, coords, edges,
                   comment=f"Cross section of the thin T-shaped domain, eps = {eps!r}, k = {k}.\n"
                           "Generated by reproducibles_python/thin_plate_eigenvalues/thin.py.")
  return len(edges)


class ThinAveraging:
  """The averaging map Pi_eps applied to discrete eigenfunctions of the cross section.

  For plate i, p_i(r) is the average of w over the thickness at distance r from the junction,
  extended to (0, 2 eps) by reflection at r = 2 eps (extension="reflection"), and r_N is the
  average over the junction. extension="zero" extends by zero instead, which is not uniformly
  bounded in H^1 (only used for comparison). The inner products of H_0 are evaluated without the
  common factor int cos^2(k pi s) ds.
  """

  def __init__(self, hdg, eps, n_quad, extension="reflection"):
    assert extension in ("reflection", "zero")
    self.eps, self.reflect = eps, extension == "reflection"
    self.abscissas, weights_1d = common.gauss(n_quad)
    self.n_quad = n_quad
    points = hdg.bulk_points(list(self.abscissas)).reshape(hdg.n_hyEdges(), n_quad, n_quad, 2)
    # points[e, b, a] with a (x-direction) running fastest
    lower, upper = points[:, 0, 0], points[:, -1, -1]
    centers = 0.5 * (lower + upper)
    self.region = np.array([region(x, y, eps) for x, y in centers])
    measures = hdg.hyEdge_measures()
    self.weights = measures[:, None, None] * np.outer(weights_1d, weights_1d)[None, :, :]
    self.points, self.hdg = points, hdg
    # radial direction: x for plates 1 and 2, y for plate 3; thickness: the other one
    self.columns = []
    for plate in range(3):
      cells = np.flatnonzero(self.region == plate)
      if plate < 2:  # cells of a column share their x-interval
        key = np.round(np.minimum(abs(lower[cells, 0]), abs(upper[cells, 0])) / eps, 10)
      else:
        key = np.round(lower[cells, 1] / eps, 10)
      columns = {}
      for cell, value in zip(cells, key):
        columns.setdefault(value, []).append(cell)
      self.columns.append([np.array(columns[value]) for value in sorted(columns)])
    self.junction = np.flatnonzero(self.region == 3)
    self.junction_area = 39. / 4. * eps * eps
    self.weights_1d = weights_1d

  def averages(self, x):
    """(r, weights, p_i(r)) per plate at the Gauss points and the junction average."""
    values = self.hdg.bulk_values(list(x), list(self.abscissas)).reshape(-1, self.n_quad,
                                                                          self.n_quad)
    result = []
    for plate in range(3):
      thickness = 2. * PLATE_HALF_THICKNESS[plate] * self.eps
      rs, ws, ps = [], [], []
      for cells in self.columns[plate]:
        pts = self.points[cells[0]]
        if plate < 2:
          r = np.abs(pts[0, :, 0])  # along a, the x-direction
          length = abs(pts[0, -1, 0] - pts[0, 0, 0]) / (self.abscissas[-1] - self.abscissas[0])
          # integrate over b (y-direction) for every a
          p = sum((values[c] * self.weights[c]).sum(axis=0) / self.weights_1d for c in cells)
          p = p / length / thickness
        else:
          r = pts[:, 0, 1]  # along b, the y-direction
          length = abs(pts[-1, 0, 1] - pts[0, 0, 1]) / (self.abscissas[-1] - self.abscissas[0])
          p = sum((values[c] * self.weights[c]).sum(axis=1) / self.weights_1d for c in cells)
          p = p / length / thickness
        rs.append(r)
        ws.append(length * self.weights_1d)
        ps.append(p)
      result.append((np.concatenate(rs), np.concatenate(ws), np.concatenate(ps)))
    junction = sum((values[c] * self.weights[c]).sum() for c in self.junction) / self.junction_area
    return result, junction

  def inner(self, averaged1, averaged2):
    """<Pi_eps w_1, Pi_eps w_2>_{H_0} (radial part) of two averaged functions."""
    value = 0.
    for plate in range(3):
      r, w, p1 = averaged1[0][plate]
      p2 = averaged2[0][plate][2]
      collar = (r < 4. * self.eps) & self.reflect  # reflected to (0, 2 eps)
      value += limit.D[plate] * (np.dot(w, p1 * p2) + np.dot(w[collar], (p1 * p2)[collar]))
    return value + limit.MU * averaged1[1] * averaged2[1]

  def inner_with_mode(self, averaged, mode):
    """<Pi_eps w, u>_{H_0} (radial part) for a limit eigenfunction u."""
    value = 0.
    for plate in range(3):
      r, w, p = averaged[0][plate]
      collar = (r < 4. * self.eps) & self.reflect
      value += limit.D[plate] * (np.dot(w, p * mode.radial(plate, r)) +
                                 np.dot(w[collar], p[collar] *
                                        mode.radial(plate, 4. * self.eps - r[collar])))
    return value + limit.MU * averaged[1] * mode.junction_value


def solve(eps, k, degree, level=1, n_eigs=3, corner_levels=4, directory=None, first=None,
          largest=1. / 16., tau=1.):
  """Eigenpairs of the cross-section problem with k^2 pi^2 reaction (k as in cos(k pi s))."""
  directory = directory or common.OUTPUT
  os.makedirs(directory, exist_ok=True)
  tag = f"{eps:.6e}_k{k}_c{corner_levels}_f{(first or eps / 2.) / eps:g}_l{largest:g}"
  filename = os.path.join(directory, f"cross_section_{tag}.geo")
  if not os.path.isfile(filename):
    write_domain(filename, eps, k, corner_levels, first, largest)
  hdg = common.problem_class(degree, 2)(filename, tau)
  hdg.set_refinement(level)
  disc = common.Discretization(hdg)
  values, vectors = disc.eigenpairs(n_eigs)
  return hdg, disc, values, vectors


# The thin-domain eigenvalues of the clusters of limit.clusters(): (k, indices in the k-problem).
CLUSTER_SECTORS = {"lambda1": (0, [0]), "lambda23": (0, [1, 2]), "lambda4": (1, [0])}


def thin_eigenpairs(eps, degree, level=1, **kwargs):
  """Eigenvalues and eigenvectors of the k = 0 and k = 1 problems."""
  return {k: solve(eps, k, degree, level, n_eigs=4 if k == 0 else 2, **kwargs) for k in (0, 1)}


def errors(eps, sectors, degree):
  """Eigenvalue errors and averaged eigenspace gaps of the clusters of limit.clusters()."""
  result = {}
  for name, cluster in limit.clusters().items():
    k, indices = CLUSTER_SECTORS[name]
    hdg, disc, values, vectors = sectors[k]
    exact = limit.eigenvalue(cluster)
    modes = cluster["modes"]
    result[name] = {
      "eigenvalues": [float(values[i]) for i in indices],
      "eigenvalue_error": float(max(abs(values[i] - exact) for i in indices)),
    }
    for extension, key in [("reflection", "gap"), ("zero", "gap_zero_extension")]:
      averaging = ThinAveraging(hdg, eps, degree + 6, extension)
      averaged = [averaging.averages(vectors[:, i]) for i in indices]
      gram_eps = np.array([[averaging.inner(a, b) for b in averaged] for a in averaged])
      cross = np.array([[averaging.inner_with_mode(a, mode) for mode in modes] for a in averaged])
      result[name][key] = common.gap(gram_eps, limit.gram(modes, limit.radial_inner), cross)
  return result


def _column_polynomials(r, p, n_quad, degree):
  """The thickness average p(r) on every column (n_quad consecutive Gauss points) as Legendre
  coefficients on its interval; p is a polynomial of degree at most `degree` on every column."""
  abscissas = common.gauss(n_quad)[0]
  pieces = []
  for first in range(0, len(r), n_quad):
    r_col, p_col = r[first:first + n_quad], p[first:first + n_quad]
    length = (r_col.max() - r_col.min()) / (abscissas[-1] - abscissas[0])
    start = r_col.min() - length * abscissas[0]
    t = 2. * (r_col - start) / length - 1.
    pieces.append((start, start + length, np.polynomial.legendre.legfit(t, p_col, degree)))
  return pieces


def _evaluate_columns(pieces, r):
  values = np.full(len(r), np.nan)
  for start, end, coeffs in pieces:
    inside = (r > start - 1e-14) & (r < end + 1e-14) & np.isnan(values)
    values[inside] = np.polynomial.legendre.legval(2. * (r[inside] - start) / (end - start) - 1.,
                                                   coeffs)
  assert not np.isnan(values).any(), "The fine mesh must refine the coarse mesh along the plates!"
  return values


def mesh_comparison(eps, coarse, fine, degree):
  """Compare the averaged eigenspaces Pi_eps M_eps of two discretizations of the cross sections.

  coarse and fine are results of thin_eigenpairs (of the same degree), where the fine mesh refines
  the coarse one. For every cluster, the result contains the gap between the two averaged
  eigenspaces (with respect to H_0) and the difference of the gaps g_eps to the limit eigenspace.
  The thickness averages of the coarse discretization are evaluated at the quadrature points of the
  fine one, which is exact since they are polynomials on the columns of the coarse mesh.
  """
  n_quad = degree + 6
  errors_coarse, errors_fine = errors(eps, coarse, degree), errors(eps, fine, degree)
  result = {}
  for name, (k, indices) in CLUSTER_SECTORS.items():
    av_coarse = ThinAveraging(coarse[k][0], eps, n_quad)
    av_fine = ThinAveraging(fine[k][0], eps, n_quad)
    f_coarse = [av_coarse.averages(coarse[k][3][:, i]) for i in indices]
    f_fine = [av_fine.averages(fine[k][3][:, i]) for i in indices]

    def cross(a, b):
      value = 0.
      for plate in range(3):
        r_a, _, p_a = a[0][plate]
        r, w, p = b[0][plate]
        q = _evaluate_columns(_column_polynomials(r_a, p_a, n_quad, degree), r)
        collar = r < 4. * eps  # reflected to (0, 2 eps)
        value += limit.D[plate] * (np.dot(w, q * p) + np.dot(w[collar], (q * p)[collar]))
      return value + limit.MU * a[1] * b[1]

    gram_coarse = np.array([[av_coarse.inner(a, b) for b in f_coarse] for a in f_coarse])
    gram_fine = np.array([[av_fine.inner(a, b) for b in f_fine] for a in f_fine])
    mixed = np.array([[cross(a, b) for b in f_fine] for a in f_coarse])
    result[name] = {
      "eigenspace_gap": common.gap(gram_coarse, gram_fine, mixed),
      "gap_difference": abs(errors_coarse[name]["gap"] - errors_fine[name]["gap"]),
      "gap_zero_extension_difference": abs(errors_coarse[name]["gap_zero_extension"] -
                                           errors_fine[name]["gap_zero_extension"]),
    }
  return result


def global_order_ok(sectors, n_check=4):
  """Whether the smallest n_check eigenvalues of the thin domain are ordered as in the limit."""
  merged = sorted([(v, k) for k in (0, 1) for v in sectors[k][2]])[:n_check]
  return [k for _, k in merged] == [0, 0, 0, 1]


if __name__ == "__main__":
  import sys
  degree = int(sys.argv[1]) if len(sys.argv) > 1 else 3
  for eps in [1. / 8., 1. / 16., 1. / 32., 1. / 64.]:
    sectors = thin_eigenpairs(eps, degree)
    res = errors(eps, sectors, degree)
    print(f"eps={eps:.4e} n={sectors[0][1].n_unknowns()} order_ok={global_order_ok(sectors)}",
          {k: (v["eigenvalues"], v["eigenvalue_error"], v["gap"]) for k, v in res.items()})
