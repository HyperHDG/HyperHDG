#!/usr/bin/env python3
"""A network of Timoshenko beams on the unit square, under a force or a displacement of a side.

    python3 examples/beam_network_nb_jit.py                        # the grid of domains/grid_7x7.geo
    python3 examples/beam_network_nb_jit.py --experiment tensile
    python3 examples/beam_network_nb_jit.py network.geo.h5 -o beam.vtkhdf

Three experiments:
    force     all four sides clamped, the force -e_z on every beam (the default)
    tensile   the left side clamped, the right side clamped and moved by 0.1 e_x
    bending   the left side clamped, the right side held in z only and moved by 0.1 e_z
--force, --displacement and --clamp change their force, the displacement of the right side and the
clamped sides (--clamp clamps all components). Clamped nodes have zero rotation and the
displacement that interpolates the displacements of the four corners of the square bilinearly
(zero on the left, --displacement on the right); components that are not clamped are free.

The domain is a .geo file or, with HDF5, a .geo.h5 file (e.g. from
github.com/HyperHDG/fiber_network), whose node types tell the sides: 1 and 2 on xmin and xmax,
3 and 6 on ymin and ymax, summed at the corners. The C++ below is the module source, in which
DEGREE is replaced; hyperhdg.load compiles it with CMake the first time a degree is asked for and
imports it. Its parameters are static data of the module, set before the system is built. The
solution is written to the file of -o, whose ending tells the format: .vtu or, with HDF5, .vtkhdf;
ParaView reads both.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def tprint(*args):
    print(f"# [{time.strftime('%H:%M:%S')}]", *args)

try:  # optional: the Cholesky factorization of CHOLMOD, else the direct solver of SciPy
    import sksparse.cholmod as spch

    def sparse_direct_solve(matrix, rhs):
        return spch.cho_factor(matrix).solve(rhs)
    tprint("using cholmod")
except ImportError:
    sparse_direct_solve = spla.spsolve
    tprint("using superlu")

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
COMPONENTS = ["ux", "uy", "uz", "rx", "ry", "rz"]  # the unknowns of a node, in this order
EXPERIMENTS = {  # force, displacement of the right side, clamped components of the clamped sides
    "force": ((0., 0., -1.), (0., 0., 0.), dict.fromkeys(SIDES, COMPONENTS)),
    "tensile": ((0., 0., 0.), (0.01, 0., 0.), {"xmin": COMPONENTS, "xmax": COMPONENTS}),
    "bending": ((0., 0., 0.), (0., 0., 0.1), {"xmin": COMPONENTS, "xmax": ["uz"]}),
}


def node_types(sides):
    """The node types on the sides: digit 0 (x) or 1 (y) of the type in base 3 names the side."""
    return [node_type for node_type in range(1, 27)
            if any(node_type // 3**axis % 3 == side
                   for axis, side in (SIDES[name] for name in sides))]


def set_problem(module, hdg, force, displacement, clamp):
    """Set the force, the displacement of the right side and the constraints of the clamped sides.

    clamp maps a side to its clamped components. The loop constrains all unknowns of the nodes
    of the clamped node types (constrain_dirichlet_nodes), so the sides are constrained one by
    one, keeping the unknowns of their clamped components (those of a node are consecutive, in
    the order of COMPONENTS). A corner has the clamped components of both of its sides.
    """
    # Free components have no prescribed value: only the clamped ones enter the load.
    displacement = [value if COMPONENTS[d] in clamp.get("xmax", COMPONENTS) else 0.
                    for d, value in enumerate(displacement)]
    corner_displacement = [[0., 0., 0.], displacement] * 2  # corners (0,0), (1,0), (0,1), (1,1)
    constraints = []
    for side, components in clamp.items():
        module.set_problem(node_types([side]), force, corner_displacement)
        hdg.constrain_dirichlet_nodes()
        unknowns = hdg.constraints()
        clamped = [COMPONENTS.index(name) for name in components]
        constraints.append(unknowns[np.isin(unknowns % len(COMPONENTS), clamped)])
    module.set_problem(node_types(clamp), force, corner_displacement)
    hdg.set_constraints(np.unique(np.concatenate(constraints)))


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
                    help="sides clamped in all components (default: those of the experiment)")
parser.add_argument("--degree", type=int, default=2, help="polynomial degree")
parser.add_argument("-o", "--output", default="output/beam_network.vtu", help="output file")
opt = parser.parse_args()

force, displacement, clamp = EXPERIMENTS[opt.experiment]
if opt.clamp:
    clamp = dict.fromkeys(opt.clamp, COMPONENTS)

tprint("compiling ...")
module = hyperhdg.load(MODULE.replace("DEGREE", str(opt.degree)))
tprint("constructing ...")
hdg = module.Beams(opt.domain)
set_problem(module, hdg, opt.force or force, opt.displacement or displacement, clamp)
tprint("preparing ...")
hdg.prepare()

# The system: the residual at zero is minus the right-hand side, assemble writes the triplets of
# the matrix, with a unit diagonal on the constrained unknowns. It is positive definite.
size = hdg.size_of_system()
rhs = np.empty(size)
hdg.residual(np.zeros(size), rhs)
m = hdg.n_triplets(True)
rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
tprint("assembling ...")
hdg.assemble(rows, cols, vals, 0., True)
matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsc()
tprint("solving ...")
solution = sparse_direct_solve(matrix, -rhs)
tprint("solved")

print(f"beams: {hdg.n_hyEdges()}")
print(f"unknowns: {size}")
# the analytic displacement of the parameters is zero: this is the norm of the solution
print(f"norm_l2: {hdg.errors(solution)[0]:.6e}")

hdg.plot_option("filePath", opt.output)
hdg.plot_option("printFileNumber", "false")
hdg.plot_solution(solution)
print(f"plot: {opt.output}")
