#!/usr/bin/env python3
"""Timoshenko beam network waves (timowave), manufactured solution wave4, from the command line.

    python python/timowave.py -deg 3 -nx 8 -nt 64 -wave4_px -1.5707963267948966

HDG in space, s-stage Gauss collocation in time: per time step every Gauss stage
representative's condensed trace system is solved with scipy's sparse LU (the matrices are
time-constant, so they are factorized once). The (degree, stages) instantiation is compiled
the first time it is asked for and cached afterwards -- see python/README.md. Stage systems of
s >= 2 are complex.

The manufactured solution lives in python/parameters/timowave.hxx: a standing wave on the
axis-aligned cross (domains/cross2.geo). -wave4_px -pi/2 makes
its Dirichlet data vanish identically (no stiff order reduction, temporal order 2s); the
default 0 drives every tip with time-dependent Dirichlet data.

The report is one 'key: value' line per quantity (YAML); errors are maxima over all time
steps, relative to the corresponding maxima of the analytic solution.
"""

import argparse
import math
import os

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import hyperhdg  # python/hyperhdg.py, next to this script

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <HyperHDG/global_loop/hyperbolic.hxx>
#include <HyperHDG/local_solver/timowave.hxx>

#include "python/parameters/timowave.hxx"

#include <complex>

using HDG = GlobalLoop::Hyperbolic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoWave<1, 3, {degree}, {quadrature}, TestTimoWave4, double, {stages}>,
  std::vector<{scalar}> >;

NB_MODULE({module}, m)
{{
  m.def("wave4_set_px", [](double px) {{ TestTimoWave4<3, double>::px = px; }},
        "spatial phase px of the wave4 standing-wave factor cos(w*(x+y+z) + px)");
  HyperHDG::bind_python<HDG, std::string, std::vector<double> >(m, "{name}");
}}
"""


def problem(degree, stages):
    """Compile the (degree, stages) instantiation; return the module and its solver class.

    Every combination gets its own module name and thereby its own build tree in the cache,
    so sweeps over several combinations can run in parallel once each one is compiled.
    """
    name = f"TimoWave_P{degree}_S{stages}"
    code = MODULE.format(degree=degree, quadrature=2 * degree, stages=stages,
                         scalar="double" if stages == 1 else "std::complex<double>",
                         module=f"timowave_p{degree}_s{stages}", name=name)
    module = hyperhdg.load(code)
    return module, module[name]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, prefix_chars="-",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-deg", type=int, default=1, help="polynomial degree")
    p.add_argument("-stages", type=int, default=0,
                   help="Gauss stages (0 = matched to degree: (deg+1)//2)")
    p.add_argument("-tau", type=float, default=1.0, help="HDG penalty parameter")
    p.add_argument("-tau_s", type=int, default=None, choices=[1, 0, -1],
                   help="penalty scaling tau ~ h^s, overrides -tau")
    p.add_argument("-nx", type=int, default=1, help="refinement level (elements per edge)")
    p.add_argument("-nt", type=int, default=1, help="number of time steps")
    p.add_argument("-T", type=float, default=1.0, help="end time")
    p.add_argument("-domain", default=os.path.join(REPO, "domains/cross2.geo"))
    p.add_argument("-test", default="wave4", choices=["wave4"])
    p.add_argument("-wave4_px", type=float, default=0.0,
                   help="spatial phase of the wave4 standing-wave factor (-pi/2: homogeneous "
                        "Dirichlet data)")
    p.add_argument("-compile_only", action="store_true",
                   help="compile the (deg, stages) instantiation and exit")
    return p.parse_args()


def main():
    opt = parse_args()
    stages = opt.stages if opt.stages != 0 else (opt.deg + 1) // 2
    module, cls = problem(opt.deg, stages)
    if opt.compile_only:
        return

    tau = opt.tau
    if opt.tau_s is not None:
        tau = float(opt.nx) ** -opt.tau_s  # h = 1/nx on the unit-length arms
    dt = opt.T / opt.nt

    module.wave4_set_px(opt.wave4_px)
    hdg = cls(opt.domain, [tau, dt])
    hdg.set_refinement(opt.nx)

    print(f"timowave_test: {opt.test}")
    print(f"wave4_px: {opt.wave4_px:g}")
    print(f"poly_deg: {opt.deg}")
    print(f"tau: {tau:.5e}")
    print(f"tau_s: {opt.tau_s if opt.tau_s is not None else 0}")
    print(f"nt: {opt.nt}")
    print(f"nx: {opt.nx}")
    print(f"dt: {dt:.5e}")
    print(f"T: {opt.T:.5e}")
    print(f"gauss_stages: {cls.n_gauss_stages()}")

    n = hdg.size_of_system()
    zero = hdg.zero_vector()

    # running maxima over all time steps; components: 0 = state L2, 1 = trace (length-weighted
    # skeleton norm), 2 = dual pair (n, m) L2
    err = np.zeros(3)
    nrm = np.zeros(3)

    def record(trace, t):
        np.maximum(err, hdg.errors(trace, t)[:3], out=err)
        np.maximum(nrm, hdg.norms(trace, t)[:3], out=nrm)

    hdg.make_initial(zero, 0.0)
    record(zero, 0.0)
    print(f"e_abs0: {err[0]:.5e}")

    # condensed stage operators are time-constant: assemble and factorize once per representative
    lu = []
    for rep in range(cls.n_gauss_reps()):
        rows, cols, vals = hdg.trace_to_flux_mat_stage(rep, 0.0)
        lu.append(spla.splu(sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()))

    trace = zero
    for step in range(1, opt.nt + 1):
        t = step * dt
        for rep, factor in enumerate(lu):
            rhs = -np.asarray(hdg.residual_flux_stage(zero, rep, t))
            trace = list(factor.solve(rhs))
            hdg.set_data_stage(trace, rep, t)  # stash the stage's local solution and trace
        hdg.finalize_step()  # combine the stages into the new state (and endpoint trace)
        record(trace, t)

    # relative errors; an analytic norm of zero (e.g. no analytic dual) keeps the absolute one
    rel = np.divide(err, nrm, out=err.copy(), where=nrm > 0)
    print(f"e_abs: {err[0]:.5e}")
    print(f"n_abs: {nrm[0]:.5e}")
    print(f"n_trace: {nrm[1]:.5e}")
    print(f"n_dual: {nrm[2]:.5e}")
    print(f"e_rel: {rel[0]:.5e}")
    print(f"e_trace: {rel[1]:.5e}")
    print(f"e_dual: {rel[2]:.5e}")


if __name__ == "__main__":
    main()
