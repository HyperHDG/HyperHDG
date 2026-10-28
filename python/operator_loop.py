#!/usr/bin/env python3
"""Solve one LDG-H diffusion problem with the local-operator loop or with GlobalLoop::Elliptic.

    python python/operator_loop.py --dim 2 --degree 3 --n 16 --loop operator --method assembled

Prints the result as YAML to stdout. solve(opt) is the whole driver: it reads the fields of
Options from any object that has them (argparse namespace, dataclass, SimpleNamespace) and
returns a Result; main() only parses the command line and dumps the result. Parameter sweeps are the shell's business
(python/operator_loop.sh with GNU parallel); python/operator_loop_convergence.py is a pure python
example; python/operator_loop_check.py compares both loops to roundoff.

The (dim, degree) instantiation of both loops is compiled on first use (python/hyperhdg.py).
"""

import argparse
import sys
import time
from dataclasses import asdict, dataclass, fields

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import yaml

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
  LocalSolver::DiffusionOperator<{dim}, {degree}, {quadrature}, TestParametersSinEllipt, double> >;

NB_MODULE({name}, m)
{{
  HyperHDG::bind_python<Elliptic, std::vector<unsigned int> >(m, "Elliptic");
  HyperHDG::bind_python<OpLoop, std::vector<unsigned int> >(m, "OpLoop");
}}
"""


@dataclass
class Options:
  dim: int = 2
  degree: int = 1
  n: int = 4                  # elements per direction of the unit cube
  tau: float = 1.0            # HDG penalty
  loop: str = "operator"      # operator | elliptic
  method: str = "assembled"   # assembled (scipy CG on the matrix) | matrixfree (CG on the product)
  prepare: bool = True        # operator: keep the operator store; elliptic: factorize_once
  rtol: float = 1e-10         # CG tolerance
  repeat: int = 1             # timings: best of this many runs


@dataclass
class Result:
  opt: dict                   # the options the result was computed with
  dofs: int
  edges: int
  operators: int              # stored operators (0 for elliptic / streaming)
  error: float                # L2 error of the primal unknown
  iterations: int             # CG iterations
  t_build: float              # construct + constraints [s]
  t_setup: float              # assembled: matrix; matrixfree: prepare [s]
  t_solve: float              # CG [s]
  t_apply: float              # one matrix-free product [s]
  t_errors: float             # error evaluation [s]


def module(dim, degree):
  name = f"hdg_oploop_d{dim}_p{degree}"
  return hyperhdg.load(MODULE.format(dim=dim, degree=degree, quadrature=2 * degree, name=name),
                       hash_name=False)


def problem(opt):
  """The problem object of the chosen loop, Dirichlet dofs constrained."""
  mod = module(opt.dim, opt.degree)
  if opt.loop == "elliptic":
    return mod.Elliptic([opt.n] * opt.dim, opt.tau)
  hdg = mod.OpLoop([opt.n] * opt.dim, opt.tau)
  hdg.constrain_dirichlet_nodes()
  return hdg


def timed(fun, repeat):
  best, result = float("inf"), None
  for _ in range(repeat):
    start = time.perf_counter()
    result = fun()
    best = min(best, time.perf_counter() - start)
  return best, result


def solve(opt):
  t_build, hdg = timed(lambda: problem(opt), opt.repeat)
  elliptic = opt.loop == "elliptic"
  n = hdg.size_of_system()
  zero = np.zeros(n)

  # the operator and its right-hand side; Elliptic masks Dirichlet rows, the prototype keeps
  # identity rows there, so both systems have the same solution. Elliptic allocates its results,
  # the prototype writes into arrays the caller owns (y, triplets).
  if elliptic:
    apply = lambda x: hdg.trace_to_flux(x)
    rhs = -hdg.residual_flux(zero)
    triplets = lambda: hdg.trace_to_flux_mat(0., opt.prepare)
  else:
    y = np.empty(n)
    def apply(x):
      hdg.apply(x, y)
      return y
    hdg.residual(zero, rhs := np.empty(n))
    rhs *= -1.
    m = hdg.n_triplets(True)
    rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
    def triplets():
      hdg.assemble(rows, cols, vals, 0., True)
      return rows, cols, vals
    if opt.prepare:
      hdg.prepare()

  if opt.method == "assembled":
    t_setup, matrix = timed(
      lambda: sp.coo_matrix(((t := triplets())[2], (t[0], t[1])), shape=(n, n)).tocsr(),
      opt.repeat)
    operator = matrix
  else:
    t_setup = timed(lambda: hdg.prepare(), opt.repeat)[0] if opt.prepare and not elliptic else 0.
    # scipy may keep the returned vector, so the matrix-free matvec hands out a copy
    operator = spla.LinearOperator((n, n), matvec=lambda x: apply(x).copy(), dtype=float)

  iterations = [0]
  def count(_): iterations[0] += 1
  t_solve, (solution, info) = timed(
    lambda: spla.cg(operator, rhs, rtol=opt.rtol, maxiter=100 * n, callback=count), 1)
  if info != 0:
    raise RuntimeError(f"CG did not converge (info = {info})")

  t_apply = timed(lambda: apply(solution), opt.repeat)[0]
  t_errors, error = timed(lambda: hdg.errors(solution)[0], opt.repeat)

  return Result(opt=dict(vars(opt)), dofs=n, edges=opt.n ** opt.dim,
                operators=0 if elliptic else hdg.n_operators(), error=float(error),
                iterations=iterations[0], t_build=t_build, t_setup=t_setup, t_solve=t_solve,
                t_apply=t_apply, t_errors=t_errors)


def parse_args(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  for f in fields(Options):
    if f.type is bool:
      parser.add_argument(f"--{f.name}", type=lambda v: v.lower() in ("1", "true", "yes"),
                          default=f.default, metavar="BOOL")
    else:
      parser.add_argument(f"--{f.name}", type=f.type, default=f.default)
  return parser.parse_args(argv)


def main():
  result = solve(parse_args())
  print(yaml.safe_dump(asdict(result), sort_keys=False), end="")


if __name__ == "__main__":
  sys.exit(main())
