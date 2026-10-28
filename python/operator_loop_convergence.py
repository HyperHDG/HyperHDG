#!/usr/bin/env python3
"""Convergence of the local-operator loop in pure python: errors and rates over n.

    python python/operator_loop_convergence.py --dim 2 --degree 3 --levels 5 [--loop elliptic]

The shell equivalent (python/operator_loop.sh) runs the same solve() in parallel via GNU parallel;
this script needs nothing but python and is the quickest way to see the method converge.
"""

import argparse
from types import SimpleNamespace

import numpy as np

import operator_loop


def main():
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--dim", type=int, default=2)
  parser.add_argument("--degree", type=int, default=1)
  parser.add_argument("--loop", default="operator")
  parser.add_argument("--levels", type=int, default=5)
  opt = parser.parse_args()

  print(f"{'n':>5} {'dofs':>8} {'error':>14} {'rate':>6}")
  previous = None
  for level in range(opt.levels):
    r = operator_loop.solve(SimpleNamespace(**{**vars(operator_loop.Options()), **vars(opt), "n": 2 ** level}))
    rate = "-" if previous is None else f"{np.log2(previous / r.error):.2f}"
    print(f"{r.opt['n']:>5} {r.dofs:>8} {r.error:>14.6e} {rate:>6}")
    previous = r.error


if __name__ == "__main__":
  main()
