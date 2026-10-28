#!/usr/bin/env python3
"""Prototype of the local-operator global loop, checked against GlobalLoop::Elliptic.

    python python/operator_loop.py check          # consistency: matrix-free, matrices, errors
    python python/operator_loop.py convergence    # error ladders of both loops side by side
    python python/operator_loop.py bench          # runtimes: assembly, matrix-free products, CG

Both loops run LDG-H diffusion (TestParametersSinEllipt) on the unit cube with Cubic topology.
Elliptic uses LocalSolver::Diffusion; the prototype uses GlobalLoop::LocalOperatorLoop with
LocalSolver::DiffusionOperator, which derives from Diffusion and only adds the operator protocol.
The C++ is compiled on first use via python/hyperhdg.py (one module per dimension and degree).

Every command takes --dim and --degree to run a single case (one process per instantiation), so
cases can run as independent processes; python/operator_loop.sh does that with GNU parallel
(compiles in parallel, timings serial). Without --dim/--degree all cases run in this process.

Options: --openmp compiles the module with OpenMP (the prototype's sweeps then run in parallel;
Elliptic stays serial), --threads N sets OMP_NUM_THREADS for the bench.
"""

import argparse
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import hyperhdg  # python/hyperhdg.py, next to this script

MODULE = """\
#include <HyperHDG/bind_python.hxx>

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/global_loop/local_operator_loop.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_operator.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

using Topo = Topology::Cubic<{dim}, {dim}>;
using Geom = Geometry::UnitCube<{dim}, {dim}, double>;
using Node = NodeDescriptor::Cubic<{dim}, {dim}>;

using Elliptic = GlobalLoop::Elliptic<
  Topo, Geom, Node,
  LocalSolver::Diffusion<{dim}, {degree}, {quadrature}, TestParametersSinEllipt, double>,
  std::vector<double> >;
using OpLoop = GlobalLoop::LocalOperatorLoop<
  Topo, Geom, Node,
  LocalSolver::DiffusionOperator<{dim}, {degree}, {quadrature}, TestParametersSinEllipt, double>,
  std::vector<double> >;

NB_MODULE({name}, m)
{{
  HyperHDG::bind_python<Elliptic, std::vector<unsigned int> >(m, "Elliptic");
  HyperHDG::bind_python<OpLoop, std::vector<unsigned int> >(m, "OpLoop");
}}
"""

OPENMP = False


def module(dim, degree):
  name = f"hdg_oploop_d{dim}_p{degree}" + ("_omp" if OPENMP else "")
  code = MODULE.format(dim=dim, degree=degree, quadrature=2 * degree, name=name)
  args = ["-DCMAKE_CXX_FLAGS=-fopenmp", "-DCMAKE_MODULE_LINKER_FLAGS=-fopenmp"] if OPENMP else []
  return hyperhdg.load(code, hash_name=False, cmake_args=args)


def problems(dim, degree, n_elements, tau=1.0):
  mod = module(dim, degree)
  ell = mod.Elliptic([n_elements] * dim, tau)
  op = mod.OpLoop([n_elements] * dim, tau)
  op.constrain_dirichlet_nodes()
  return ell, op


def csr(coo, n):
  rows, cols, vals = coo
  return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def solve_cg(matrix, rhs, rtol=1e-12):
  solution, info = spla.cg(matrix, rhs, rtol=rtol, maxiter=100000)
  assert info == 0, f"CG did not converge (info = {info})"
  return solution


def timed(fun, repeat=3):
  best = float("inf")
  for _ in range(repeat):
    start = time.perf_counter()
    result = fun()
    best = min(best, time.perf_counter() - start)
  return best, result


# ---------------------------------------------------------------------------------------------
def check(dim, degree, n_elements, out):
  ell, op = problems(dim, degree, n_elements)
  n = ell.size_of_system()
  rng = np.random.default_rng(1)
  x = rng.standard_normal(n)
  constrained = np.asarray(op.constraints(), dtype=int)
  free = np.ones(n, dtype=bool)
  free[constrained] = False

  def report(label, value, tol):
    status = "ok" if value <= tol else "FAIL"
    line = f"  {label:<52s} {value:10.3e}  {status}"
    print(line)
    out.write(line + "\n")
    return value <= tol

  passed = True
  # matrix-free product: Elliptic zeroes Dirichlet rows, the prototype puts x_i there
  y_ell = np.asarray(ell.trace_to_flux(list(x)))
  y_op = np.asarray(op.apply(list(x)))
  scale = np.abs(y_ell).max()
  passed &= report("apply vs trace_to_flux (free rows, rel)",
                   np.abs(y_ell - y_op)[free].max() / scale, 1e-13)
  passed &= report("apply on constrained rows = x_i",
                   np.abs(y_op[constrained] - x[constrained]).max(), 0.)

  r_ell = np.asarray(ell.residual_flux(list(x)))
  r_op = np.asarray(op.residual(list(x)))
  passed &= report("residual vs residual_flux (free rows, rel)",
                   np.abs(r_ell - r_op)[free].max() / np.abs(r_ell).max(), 1e-13)

  a_ell = csr(ell.trace_to_flux_mat(), n)
  a_op = csr(op.assemble(), n)
  passed &= report("assemble vs trace_to_flux_mat (rel)",
                   abs(a_ell - a_op).max() / abs(a_ell).max(), 1e-13)

  # prepared store: identical results, bit for bit
  op.prepare()
  y_prep = np.asarray(op.apply(list(x)))
  a_prep = csr(op.assemble(), n)
  # (shared operators come from a representative hyperedge: equal up to roundoff)
  passed &= report("prepared (shared) vs streaming apply (rel)",
                   np.abs(y_prep - y_op).max() / scale, 1e-13)
  passed &= report("prepared (shared) vs streaming assemble (rel)",
                   abs(a_prep - a_op).max() / abs(a_ell).max(), 1e-13)
  # the assembled matrix is apply on unit vectors, by construction
  # (dgetrs with n_loc right-hand sides vs one: same factors, last bits may differ)
  e = np.zeros(n)
  e[n // 3] = 1.
  col = np.asarray(op.apply(list(e)))
  passed &= report("assemble column = apply(unit vector) (free rows, rel)",
                   np.abs(a_prep[:, n // 3].toarray().ravel() - col)[free].max() / scale, 1e-14)

  # solve both ways, compare solutions and errors
  rhs_ell = -np.asarray(ell.residual_flux(list(np.zeros(n))))
  sol_ell = solve_cg(a_ell, rhs_ell)
  a_id = csr(op.assemble(0., True), n)
  rhs_op = -np.asarray(op.residual(list(np.zeros(n))))
  sol_op = solve_cg(a_id, rhs_op)
  passed &= report("CG solution, identity rows vs Elliptic (rel)",
                   np.abs(sol_ell - sol_op).max() / np.abs(sol_ell).max(), 1e-9)
  err_ell, err_op = ell.errors(list(sol_ell))[0], op.errors(list(sol_op))[0]
  passed &= report("L2 error Elliptic vs prototype (rel)",
                   abs(err_ell - err_op) / err_ell, 1e-9)
  line = f"  operators stored: {op.n_operators()} for {n} dofs, error {err_op:.6e}"
  print(line)
  out.write(line + "\n")
  op.release()
  return passed


def run_check(opt):
  os.makedirs("output", exist_ok=True)
  all_passed = True
  with open("output/operator_loop_check.txt", "w") as out:
    for dim, degree, n in [(1, 1, 8), (1, 3, 8), (2, 1, 4), (2, 2, 4), (2, 3, 3), (3, 1, 2),
                           (3, 2, 2)]:
      head = f"dim {dim}, degree {degree}, {n}^{dim} elements"
      print(head)
      out.write(head + "\n")
      all_passed &= check(dim, degree, n, out)
    verdict = "ALL CHECKS PASSED" if all_passed else "SOME CHECKS FAILED"
    print(verdict)
    out.write(verdict + "\n")
  return 0 if all_passed else 1


# ---------------------------------------------------------------------------------------------
CONVERGENCE_HEADER = f"{'dim':>4} {'deg':>4} {'n':>5} {'dofs':>8} {'err Elliptic':>14} " \
                     f"{'err prototype':>14} {'rel diff':>10} {'rate':>6}"


def convergence_case(dim, degree, out):
  previous = None
  for it in range(6 if dim == 1 else 5):
    n = 2 ** it
    ell, op = problems(dim, degree, n)
    size = ell.size_of_system()
    sol_ell = solve_cg(csr(ell.trace_to_flux_mat(), size),
                       -np.asarray(ell.residual_flux(ell.zero_vector())))
    op.prepare()
    sol_op = solve_cg(csr(op.assemble(0., True), size),
                      -np.asarray(op.residual(op.zero_vector())))
    err_ell, err_op = ell.errors(list(sol_ell))[0], op.errors(list(sol_op))[0]
    rate = "-" if previous is None else f"{np.log2(previous / err_op):.2f}"
    line = f"{dim:>4} {degree:>4} {n:>5} {size:>8} {err_ell:>14.6e} {err_op:>14.6e} " \
           f"{abs(err_ell - err_op) / err_ell:>10.1e} {rate:>6}"
    print(line)
    out.write(line + "\n")
    previous = err_op


def run_convergence(opt):
  os.makedirs("output", exist_ok=True)
  cases = [(opt.dim, opt.degree)] if opt.dim else [(d, p) for d in (1, 2) for p in (1, 2, 3)]
  suffix = f"_d{opt.dim}_p{opt.degree}" if opt.dim else ""
  with open(f"output/operator_loop_convergence{suffix}.txt", "w") as out:
    print(CONVERGENCE_HEADER)
    out.write(CONVERGENCE_HEADER + "\n")
    for dim, degree in cases:
      convergence_case(dim, degree, out)
      print()
      out.write("\n")
  return 0


# ---------------------------------------------------------------------------------------------
def bench_case(dim, degree, n_elements, out):
  ell, op = problems(dim, degree, n_elements)
  n = ell.size_of_system()
  n_edges = n_elements ** dim
  head = f"dim {dim}, degree {degree}, {n_edges} hyperedges, {n} dofs"
  print(head)
  out.write(head + "\n")
  rows = []

  def row(label, seconds, extra=""):
    line = f"  {label:<46s} {seconds * 1e3:10.1f} ms  {extra}"
    print(line)
    out.write(line + "\n")
    rows.append((label, seconds))

  t, _ = timed(lambda: ell.trace_to_flux_mat())
  row("Elliptic trace_to_flux_mat (per edge)", t)
  t, _ = timed(lambda: ell.trace_to_flux_mat(0., True))
  row("Elliptic trace_to_flux_mat (factorize_once)", t)
  op.release()
  t, _ = timed(lambda: op.assemble())
  row("prototype assemble, streaming", t)
  t, _ = timed(lambda: op.prepare())
  row("prototype prepare", t)
  t, _ = timed(lambda: op.assemble())
  row("prototype assemble, prepared", t)

  x = list(np.random.default_rng(0).standard_normal(n))
  reps = 10
  t, _ = timed(lambda: [ell.trace_to_flux(x) for _ in range(reps)])
  row(f"Elliptic trace_to_flux x{reps}", t)
  op.release()
  t, _ = timed(lambda: [op.apply(x) for _ in range(reps)])
  row(f"prototype apply x{reps}, streaming", t)
  op.prepare()
  t, _ = timed(lambda: [op.apply(x) for _ in range(reps)])
  row(f"prototype apply x{reps}, prepared", t)
  t, _ = timed(lambda: ell.errors(x))
  row("Elliptic errors", t)
  t, _ = timed(lambda: op.errors(x))
  row("prototype errors (prepared)", t)
  op.release()
  print()
  out.write("\n")


BENCH_CASES = [(2, 1, 64), (2, 3, 32), (3, 1, 8), (3, 2, 6), (3, 3, 4)]


def run_bench(opt):
  os.makedirs("output", exist_ok=True)
  cases = [c for c in BENCH_CASES if not opt.dim or (c[0], c[1]) == (opt.dim, opt.degree)]
  suffix = (f"_d{opt.dim}_p{opt.degree}" if opt.dim else "") + \
           (f"_omp{opt.threads}" if OPENMP else "")
  with open(f"output/operator_loop_bench{suffix}.txt", "w") as out:
    note = f"OpenMP: {OPENMP}, threads: {opt.threads if OPENMP else 1}"
    print(note)
    out.write(note + "\n")
    for dim, degree, n in cases:
      bench_case(dim, degree, n, out)
  return 0


def run_compile(opt):
  module(opt.dim, opt.degree)
  return 0


def main():
  global OPENMP
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("command", choices=["check", "convergence", "bench", "compile"])
  parser.add_argument("--dim", type=int, help="run one case: dimension")
  parser.add_argument("--degree", type=int, help="run one case: polynomial degree")
  parser.add_argument("--openmp", action="store_true", help="compile the module with OpenMP")
  parser.add_argument("--threads", type=int,
                      default=int(os.environ.get("OMP_NUM_THREADS", os.cpu_count())))
  opt = parser.parse_args()
  OPENMP = opt.openmp
  if OPENMP:
    os.environ["OMP_NUM_THREADS"] = str(opt.threads)
  return {"check": run_check, "convergence": run_convergence, "bench": run_bench,
          "compile": run_compile}[opt.command](opt)


if __name__ == "__main__":
  sys.exit(main())
