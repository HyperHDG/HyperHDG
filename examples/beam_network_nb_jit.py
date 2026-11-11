#!/usr/bin/env python3
"""A network of Timoshenko beams on the unit square, under a force or a displacement of a side.

    python3 examples/beam_network_nb_jit.py                        # the grid of domains/grid_7x7.geo
    python3 examples/beam_network_nb_jit.py --experiment tensile
    python3 examples/beam_network_nb_jit.py network.geo.h5 --solver cg --plot vtkhdf

Three experiments:
    force     all four sides clamped, the force -e_z on every beam (the default)
    tensile   the left side clamped, the right side clamped and moved by 0.1 e_x
    bending   the left side clamped, the right side clamped and moved by 0.1 e_z
--force, --displacement and --clamp change their force, the displacement of the right side and the
clamped sides. Clamped nodes have zero rotation and the displacement that interpolates the
displacements of the four corners of the square bilinearly (zero on the left, --displacement on
the right).

The domain is a .geo file or, with HDF5, a .geo.h5 file (e.g. from
github.com/HyperHDG/fiber_network), whose node types tell the sides: 1 and 2 on xmin and xmax,
3 and 6 on ymin and ymax, summed at the corners. The C++ below is the module source, in which
DEGREE is replaced; hyperhdg.load compiles it with CMake the first time a degree is asked for and
imports it. Its parameters are static data of the module, set before the system is built. The
solution is written to output/beam_network.vtu or, with --plot vtkhdf, to
output/beam_network.vtkhdf, which needs HDF5; ParaView reads both.
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

#include <array>
#include <vector>

// The problem is LocalSolver::TimoschenkoBeamParametersUnitSquare, whose static data is set at run
// time (set_problem): the node types that are clamped, the force on every beam and the
// displacements of the corners (0,0), (1,0), (0,1), (1,1) of the unit square.
template <unsigned int space_dimT, typename param_float_t = double>
using Problem = LocalSolver::TimoschenkoBeamParametersUnitSquare<space_dimT, param_float_t>;

using Beams = GlobalLoop::Generic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoBeam<1, 3, DEGREE, 2 * DEGREE, Problem> >;

NB_MODULE(beam_network_jit, m)
{
  HyperHDG::bind_python<Beams, std::string>(m, "Beams");
  m.def("set_problem",
        [](const std::vector<unsigned int>& dirichlet_nodes, const std::vector<double>& force,
           const std::vector<std::vector<double> >& corner_displacement)
        {
          Problem<3>::dirichlet_nodes = dirichlet_nodes;
          for (unsigned int d = 0; d < 3; ++d)
          {
            Problem<3>::force[d] = force.at(d);
            for (unsigned int corner = 0; corner < 4; ++corner)
              Problem<3>::corner_displacement[corner][d] = corner_displacement.at(corner).at(d);
          }
        });
}
"""

SIDES = {"xmin": (0, 1), "xmax": (0, 2), "ymin": (1, 1), "ymax": (1, 2)}
EXPERIMENTS = {  # force, displacement of the right side, clamped sides
    "force": ((0., 0., -1.), (0., 0., 0.), ("xmin", "xmax", "ymin", "ymax")),
    "tensile": ((0., 0., 0.), (0.1, 0., 0.), ("xmin", "xmax")),
    "bending": ((0., 0., 0.), (0., 0., 0.1), ("xmin", "xmax")),
}

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("domain", nargs="?", default="domains/grid_7x7.geo",
                    help=".geo or .geo.h5 file of the network")
parser.add_argument("--experiment", choices=list(EXPERIMENTS), default="force")
parser.add_argument("--force", type=float, nargs=3, metavar=("FX", "FY", "FZ"),
                    help="force on every beam (default: that of the experiment)")
parser.add_argument("--displacement", type=float, nargs=3, metavar=("UX", "UY", "UZ"),
                    help="displacement of the right side (default: that of the experiment)")
parser.add_argument("--clamp", nargs="+", choices=list(SIDES),
                    help="clamped sides (default: those of the experiment)")
parser.add_argument("--degree", type=int, default=2, help="polynomial degree")
parser.add_argument("--solver", choices=["direct", "cg"], default="direct",
                    help="sparse direct solve of the assembled matrix, or matrix-free CG")
parser.add_argument("--plot", choices=["vtu", "vtkhdf"], default="vtu", help="format of the plot")
opt = parser.parse_args()

force, displacement, clamp = EXPERIMENTS[opt.experiment]
force = opt.force or force
displacement = opt.displacement or displacement
clamp = opt.clamp or clamp
# The node types on the clamped sides: digit 0 (x) or 1 (y) of the type in base 3 names the side.
dirichlet_nodes = [node_type for node_type in range(1, 27)
                   if any(node_type // 3**axis % 3 == side
                          for axis, side in (SIDES[name] for name in clamp))]
zero = [0., 0., 0.]

module = hyperhdg.load(MODULE.replace("DEGREE", str(opt.degree)))
module.set_problem(dirichlet_nodes, force, [zero, displacement, zero, displacement])
hdg = module.Beams(opt.domain)
hdg.constrain_dirichlet_nodes()
hdg.prepare()

# The system: the residual at zero is minus the right-hand side. The condensed matrix is positive
# definite.
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
    def apply(x):
        y = np.empty(size)
        hdg.apply(x, y)
        return y
    solution, info = spla.cg(spla.LinearOperator((size, size), matvec=apply), -rhs,
                             rtol=1e-10, maxiter=20 * size)
    if info != 0:
        sys.exit(f"CG did not converge (info = {info})")

print(f"beams: {hdg.n_hyEdges()}, unknowns: {size}")
# the analytic displacement of the parameters is zero: this is the norm of the solution
print(f"L2 norm of the displacement: {hdg.errors(solution)[0]:.6e}")

hdg.plot_option("fileName", "beam_network")
hdg.plot_option("fileEnding", opt.plot)
hdg.plot_option("printFileNumber", "false")
hdg.plot_solution(solution)
print(f"plot: output/beam_network.{opt.plot}")
