#!/usr/bin/env python3
"""Diffusion on a hypergraph read from file, all of it from the command line.

    python python/diffusion.py domains/cross.geo --degree 2 --refine 8

The explicit-binding counterpart of reproducibles_python/fiber_network_diffusion.py: the
LDG-H system of a network of one-dimensional hyperedges, assembled by HyperHDG and solved by
scipy's CG. The C++ below is compiled the first time a (space dimension, degree) combination
is asked for and cached afterwards, so there is no build step -- see python/README.md.
"""

import argparse
import os
import re

import numpy as np
import scipy.sparse as sp

import hyperhdg  # python/hyperhdg.py, next to this script

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

using HDG = GlobalLoop::Elliptic<
  Topology::File<1, {dim}>,
  Geometry::File<1, {dim}>,
  NodeDescriptor::File<1, {dim}>,
  LocalSolver::Diffusion<1, {degree}, {quadrature}, TestParametersSinEllipt, double>,
  std::vector<double> >;

NB_MODULE(hdg_diffusion, m) {{ HyperHDG::bind_python<HDG>(m, "{name}"); }}
"""


def space_dim(domain):
    """Space dimension declared in the header of a .geo file."""
    with open(domain) as file:
        return int(re.search(r"Space_Dim\s*=\s*(\d+)", file.read()).group(1))


def problem(domain, degree, tau):
    """Compile the (space dimension, degree) instantiation and construct the problem."""
    dim = space_dim(domain)
    name = f"Diffusion_D{dim}_P{degree}"
    code = MODULE.format(dim=dim, degree=degree, quadrature=2 * degree, name=name)
    module = hyperhdg.load(code)
    return module[name](domain, tau)  # the class name is only known at runtime


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("domain", help="hypergraph file (.geo)")
    parser.add_argument("--degree", type=int, default=1, help="polynomial degree")
    parser.add_argument("--refine", type=int, default=1, help="uniform refinement level")
    parser.add_argument("--tau", type=float, default=1.0, help="HDG penalty parameter")
    opt = parser.parse_args()

    hdg = problem(opt.domain, opt.degree, opt.tau)
    hdg.set_refinement(opt.refine)

    size = hdg.size_of_system()
    rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
    rows, cols, vals = hdg.trace_to_flux_mat()
    matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsr()

    solution, info = sp.linalg.cg(matrix, rhs, rtol=1e-10)
    if info != 0:
        raise RuntimeError(f"CG did not converge (info = {info})")

    print(f"unknowns: {size}")
    print(f"error: {hdg.errors(list(solution))[0]:.5e}")

    os.makedirs("output", exist_ok=True)
    hdg.plot_option("fileName", f"{os.path.basename(opt.domain).split('.')[0]}-P{opt.degree}")
    hdg.plot_option("printFileNumber", "false")
    hdg.plot_solution(list(solution))


if __name__ == "__main__":
    main()
