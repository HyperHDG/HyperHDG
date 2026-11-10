#!/usr/bin/env python3
"""A network of Timoshenko beams in space, clamped at its Dirichlet nodes, under a uniform force.

    python3 examples/beam_network_nb_jit.py                        # the grid of domains/grid_7x7.geo
    python3 examples/beam_network_nb_jit.py network.geo.h5 --solver cg --plot vtkhdf

The domain is a .geo file or, with HDF5, a .geo.h5 file (e.g. from
github.com/HyperHDG/fiber_network). The force is -e_z on every beam, displacements and rotations
are zero at the Dirichlet nodes (LocalSolver::TimoschenkoBeamParametersDefault). The C++ below is
the module source with the degree filled in; hyperhdg.load compiles it with CMake the first time a
degree is asked for and imports it. The solution is written to output/beam_network.vtu or, with
--plot vtkhdf, to output/beam_network.vtkhdf, which needs HDF5; ParaView reads both.
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

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

using Beams = GlobalLoop::Generic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoBeam<1, 3, {degree}, {quadrature}> >;

NB_MODULE(beam_network_jit, m)
{{
  HyperHDG::bind_python<Beams, std::string>(m, "Beams");
}}
"""

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("domain", nargs="?", default="domains/grid_7x7.geo",
                    help=".geo or .geo.h5 file of the network")
parser.add_argument("--degree", type=int, default=2, help="polynomial degree")
parser.add_argument("--solver", choices=["direct", "cg"], default="direct",
                    help="sparse direct solve of the assembled matrix, or matrix-free CG")
parser.add_argument("--plot", choices=["vtu", "vtkhdf"], default="vtu", help="format of the plot")
opt = parser.parse_args()

module = hyperhdg.load(MODULE.format(degree=opt.degree, quadrature=2 * opt.degree))
hdg = module.Beams(opt.domain)
hdg.constrain_dirichlet_nodes()
hdg.prepare()

# The system: the residual at zero is minus the right-hand side. The condensed matrix is negative
# definite, so the system is solved with its negative.
size = hdg.size_of_system()
rhs = np.empty(size)
hdg.residual(np.zeros(size), rhs)
if opt.solver == "direct":
    # assemble writes the triplets of the matrix, with a unit diagonal on the Dirichlet dofs
    m = hdg.n_triplets(True)
    rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
    hdg.assemble(rows, cols, vals, 0., True)
    matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsc()
    solution = spla.spsolve(matrix, -rhs)
else:
    def minus_apply(x):
        y = np.empty(size)
        hdg.apply(x, y)
        return -y
    solution, info = spla.cg(spla.LinearOperator((size, size), matvec=minus_apply), rhs,
                             rtol=1e-10, maxiter=20 * size)
    if info != 0:
        sys.exit(f"CG did not converge (info = {info})")

print(f"beams: {hdg.n_hyEdges()}, unknowns: {size}")
# the analytic displacement of the default parameters is zero: this is the norm of the solution
print(f"L2 norm of the displacement: {hdg.errors(solution)[0]:.6e}")

hdg.plot_option("fileName", "beam_network")
hdg.plot_option("fileEnding", opt.plot)
hdg.plot_option("printFileNumber", "false")
hdg.plot_solution(solution)
print(f"plot: output/beam_network.{opt.plot}")
