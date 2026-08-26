#!/usr/bin/env python3
"""Cross-check the explicit nanobind bindings against the legacy cython path.

Both paths drive the very same C++ code (GlobalLoop::Elliptic + LocalSolver::Diffusion on a
file based hypergraph), so every quantity that crosses the language boundary must agree to
the last bit: the right hand side, the condensed system matrix, the CG solution and the
error. This script runs a small matrix of (domain, polynomial degree, refinement, tau) and
reports the maximum deviation.

    python python/check_vs_cython.py

Requires cython (the legacy path compiles at import time, ~4 s per configuration) and the
legacy build layout, i.e. a build directory named build/ holding cmake_cython.cfg.
"""

import os
import sys

import numpy as np
import scipy.sparse as sp

import diffusion  # the example under test, next to this script

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CASES = [
    # domain,                 space dim, poly degree, refinement, tau
    ("domains/simplex_1_2.geo", 2, 1, 1, 1.0),
    ("domains/simplex_1_2.geo", 2, 2, 1, 1.0),
    ("domains/simplex_1_2.geo", 2, 3, 2, 1.0),
    ("domains/injection_test.geo", 2, 2, 4, 2.0),
    ("domains/simplex_1_3.geo", 3, 1, 1, 1.0),
    ("domains/simplex_1_3.geo", 3, 3, 3, 0.5),
    ("domains/cross.geo", 3, 2, 8, 1.0),
]


def quantities(hdg, coo):
    """The numbers a python driver gets out of a global loop, for one problem.

    coo() returns the condensed system as (rows, cols, values); the two bindings spell that
    differently (see cython/elliptic.pyx for the legacy order).
    """
    size = hdg.size_of_system()
    rhs = -np.asarray(hdg.residual_flux(hdg.zero_vector()))
    rows, cols, vals = coo()
    matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsr()
    solution, info = sp.linalg.cg(matrix, rhs, rtol=1e-12)
    assert info == 0, f"CG did not converge (info = {info})"
    return size, rhs, matrix, solution, hdg.errors(list(solution))[0]


def nanobind_run(domain, dim, degree, refinement, tau):
    hdg = diffusion.problem(domain, degree, tau)
    hdg.set_refinement(refinement)
    return quantities(hdg, hdg.trace_to_flux_mat)


def cython_run(domain, dim, degree, refinement, tau):
    sys.path.append(os.path.join(REPO, "import"))
    import HyperHDG

    point = f"File<1,{dim},std::vector,Point<{dim},double> >"
    const = HyperHDG.config()
    const.global_loop = "Elliptic"
    const.topology = const.geometry = const.node_descriptor = point
    const.local_solver = \
        f"Diffusion<1,{degree},{2 * degree},TestParametersSinEllipt,double>"
    const.cython_replacements = ["string", "string"]
    const.include_files = ["reproducibles_python/parameters/diffusion.hxx"]
    hdg = HyperHDG.include(const)(domain, lsol_constr=tau)
    hdg.refine(refinement)

    def coo():  # legacy order is (cols, rows, values)
        cols, rows, vals = hdg.sparse_stiff_mat()
        return rows, cols, vals

    return quantities(hdg, coo)


def main():
    worst = 0.
    for domain, dim, degree, refinement, tau in CASES:
        os.chdir(REPO)
        new = nanobind_run(domain, dim, degree, refinement, tau)
        old = cython_run(domain, dim, degree, refinement, tau)
        assert new[0] == old[0], f"system size {new[0]} != {old[0]}"
        deviation = max(np.max(np.abs(new[1] - old[1])),
                        abs(new[2] - old[2]).max() if new[0] else 0.,
                        np.max(np.abs(new[3] - old[3])),
                        abs(new[4] - old[4]))
        worst = max(worst, deviation)
        print(f"{os.path.basename(domain):<22} P{degree} ref {refinement} tau {tau:<4} "
              f"n {new[0]:>6}  error {new[4]:.12e}  max deviation {deviation:.3e}")

    print(f"\nworst deviation over {len(CASES)} configurations: {worst:.3e}")
    return 0 if worst == 0. else 1


if __name__ == "__main__":
    sys.exit(main())
