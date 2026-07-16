#!/usr/bin/env python3
"""Minimal python port of experiments/timowave.cxx: the serial wave4 path only.

Consumes the explicitly instantiated nanobind module timowave_py (built by
python/CMakeLists.txt) and solves the condensed stage systems with scipy's sparse LU
instead of PETSc. Option names and the YAML-style report mirror the PETSc driver so the
existing output pipelines parse both.
"""

import argparse
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# compiled (poly_deg, stages) ladder; extend python/timowave_py.cxx to grow it
INSTANTIATIONS = {
    (1, 1): "TimoWave4_P1S1",
    (2, 2): "TimoWave4_P2S2",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, prefix_chars="-")
    p.add_argument("-deg", type=int, default=1, help="polynomial degree")
    p.add_argument("-stages", type=int, default=0,
                   help="Gauss stages (0 = matched to degree: (deg+1)//2)")
    p.add_argument("-tau", type=float, default=1.0, help="HDG penalty parameter")
    p.add_argument("-tau_s", type=int, default=None,
                   help="penalty exponent s with tau ~ h^s (1, 0, -1)")
    p.add_argument("-nx", type=int, default=1, help="number of refinements")
    p.add_argument("-nt", type=int, default=1, help="number of timesteps")
    p.add_argument("-T", type=float, default=1.0, help="end time")
    p.add_argument("-domain", default=os.path.join(REPO, "domains/single1.geo"))
    p.add_argument("-test", default="wave4", choices=["wave4"])
    p.add_argument("-wave4_px", type=float, default=0.0,
                   help="spatial phase of the wave4 standing-wave factor")
    p.add_argument("-print_timestep", action="store_true")
    p.add_argument("-build", default="openblas",
                   help="build preset whose python/ dir holds timowave_py")
    return p.parse_args()


def main():
    opt = parse_args()
    sys.path.insert(0, os.path.join(REPO, "build", opt.build, "python"))
    import timowave_py as hy

    stages = opt.stages if opt.stages != 0 else (opt.deg + 1) // 2
    try:
        cls = getattr(hy, INSTANTIATIONS[(opt.deg, stages)])
    except KeyError:
        sys.exit(f"no compiled instantiation for poly_deg = {opt.deg}, stages = {stages} "
                 f"(available: {sorted(INSTANTIATIONS)}; add it to python/timowave_py.cxx)")

    tau = opt.tau
    if opt.tau_s is not None:
        tau = {1: 1.0 / opt.nx, 0: 1.0, -1: float(opt.nx)}[opt.tau_s]
    dt = opt.T / opt.nt

    hy.wave4_set_px(opt.wave4_px)
    print(f"wave4_px: {opt.wave4_px:g}")

    hdg = cls(opt.domain, [tau, dt, 0.0])
    if opt.nx != 1:
        hdg.set_refinement(opt.nx)

    print(f"timowave_test: {opt.test}")
    print(f"poly_deg: {opt.deg}")
    print(f"tau: {tau:.5e}")
    print(f"tau_s: {opt.tau_s if opt.tau_s is not None else 0}")
    print(f"nt: {opt.nt}")
    print(f"nx: {opt.nx}")
    print(f"dt: {dt:.5e}")
    print(f"T: {opt.T:.5e}")
    print(f"gauss_stages: {cls.n_gauss_stages()}")

    n = hdg.n_local_dofs()  # serial: == size_of_system()
    zero = hdg.zero_vector()

    # error/norm history: index 0 = state L2, index 1 = trace norm (as in the PETSc driver)
    e_abs = n_abs = e_trace = n_trace = 0.0

    def record(span, t):
        nonlocal e_abs, n_abs, e_trace, n_trace
        e = hdg.errors(span, t)
        nn = hdg.norms(span, t)
        e_abs, n_abs = max(e[0], e_abs), max(nn[0], n_abs)
        e_trace, n_trace = max(e[1], e_trace), max(nn[1], n_trace)
        return e[0]

    hdg.make_initial([0.0] * n)
    print(f"e_abs0: {record(zero, 0.0):.5e}")

    # condensed stage operators, assembled once (time-constant), one per Gauss representative
    reps = cls.n_gauss_reps()
    lu = []
    for rep in range(reps):
        rows, cols, vals = hdg.trace_to_flux_mat_stage(rep, 0.0)
        A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
        lu.append(spla.splu(A))

    span = zero
    for i in range(1, opt.nt + 1):
        if opt.print_timestep:
            print(f"#------------ TIMESTEP {i} -------")
        t = i * dt
        for rep in range(reps):
            rhs = -np.asarray(hdg.residual_flux_stage(zero, rep, t))
            span = list(lu[rep].solve(rhs))
            hdg.set_data_stage(span, rep, t)  # stage trace zeta -> stash stage locals + trace
        hdg.finalize_step()
        record(span, t)

    e_rel = e_abs / n_abs
    # relative error in the length-weighted skeleton norm; without an analytic trace norm
    # (n_trace == 0) keep the absolute value
    if n_trace > 0:
        e_trace /= n_trace
    print(f"e_abs: {e_abs:.5e}")
    print(f"n_abs: {n_abs:.5e}")
    print(f"n_trace: {n_trace:.5e}")
    print(f"e_rel: {e_rel:.5e}")
    print(f"e_trace: {e_trace:.5e}")


if __name__ == "__main__":
    main()
