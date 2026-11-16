#!/usr/bin/env python3
"""LDG-H diffusion on the unit cube, compiled just in time for the dimension and degree asked for.

    python3 examples/diffusion_elliptic_nb_jit.py --dim 3 --degree 2

The C++ below is the module source with the dimension and degree filled in. hyperhdg.load compiles
it with CMake the first time a combination is asked for and imports it; later runs find it in
./.hyperhdg-cache. The solve is the one of examples/diffusion_elliptic_nb_precompiled.py.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import hyperhdg  # noqa: E402  (python/hyperhdg.py)

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <examples/parameters/diffusion.hxx>

using Diffusion = GlobalLoop::Generic<
  Topology::Cubic<{dim}, {dim}>,
  Geometry::UnitCube<{dim}, {dim}, double>,
  NodeDescriptor::Cubic<{dim}, {dim}>,
  LocalSolver::Diffusion<{dim}, {degree}, {quadrature}, HG<{dim}>::DiffusionElliptic, double> >;

NB_MODULE(diffusion_elliptic_jit, m)
{{
  HyperHDG::bind_python<Diffusion, std::vector<unsigned int> >(m, "Diffusion");
}}
"""

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--dim", type=int, default=2, help="dimension of the unit cube")
parser.add_argument("--degree", type=int, default=1, help="polynomial degree")
parser.add_argument("--n", type=int, default=8, help="elements per direction")
opt = parser.parse_args()

module = hyperhdg.load(MODULE.format(dim=opt.dim, degree=opt.degree, quadrature=2 * opt.degree))
hdg = module.Diffusion([opt.n] * opt.dim)
hdg.constrain_dirichlet_nodes()
hdg.prepare()

size = hdg.size_of_system()
rhs = np.empty(size)
hdg.residual(np.zeros(size), rhs)
rhs *= -1.
m = hdg.n_triplets(True)
rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
hdg.assemble(rows, cols, vals, 0., True)
matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsc()

solution = spla.spsolve(matrix, rhs)
print(f"unknowns: {size}")
print(f"error: {hdg.errors(solution)[0]:.6e}")
