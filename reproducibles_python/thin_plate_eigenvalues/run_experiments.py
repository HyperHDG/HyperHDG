#!/usr/bin/env python3
"""Numerical experiments of the manuscript on eigenvalue problems of thin plate networks.

    python reproducibles_python/thin_plate_eigenvalues/run_experiments.py [--quick] [--only ...]

Three studies (cf. README.md in this directory):
  hdg       HDG approximation of the limit problem on the hypergraph of three plates: eigenvalue
            errors and eigenspace gaps for the polynomial degrees 0, ..., 3 under refinement.
  thin      Eigenvalues of the thin domains for eps = 2^-3, ..., 2^-10 (by separation of variables
            and HDG on the cross sections): errors with respect to the limit eigenvalues and gaps
            between the averaged thin eigenspaces and the limit eigenspaces.
  combined  Differences between the thin-domain eigenvalues and the HDG eigenvalues (degree 1)
            for all combinations of eps and h (uses the results of hdg and thin).
The results are written to results/*.json, and the tables of the manuscript to results/*.tex.
"""

import argparse
import json
import math
import os
import platform
import sys
import time

import numpy as np

import common
import limit
import plates
import thin

RESULTS = os.path.join(common.HERE, "results")
CLUSTERS = ["lambda1", "lambda23", "lambda4"]

# The reference discretization of the cross sections and a finer one for the accuracy estimate.
THIN_REFERENCE = {"degree": 5, "level": 1, "corner_levels": 3, "largest": 1. / 32.}
THIN_CHECK = {"degree": 5, "level": 2, "corner_levels": 4, "largest": 1. / 32.}


def hdg_study(quick):
  # n = 2^level squares per plate side (for p = 0 and n = 2, lambda_2 and lambda_3 are not resolved)
  levels = {0: range(2, 8), 1: range(1, 7), 2: range(1, 6), 3: range(1, 5)}
  if quick:
    levels = {1: range(1, 4), 2: range(1, 3)}
  runs = []
  for degree, exponents in levels.items():
    for exponent in exponents:
      n = 2 ** exponent
      start = time.time()
      hdg, disc, values, vectors = plates.run(degree, n, n_eigs=6)
      errors = plates.errors(hdg, disc, values, vectors, degree)
      runs.append({"degree": degree, "n": n, "h": math.sqrt(2.) / n,
                   "unknowns": disc.n_unknowns(), "eigenvalues": values.tolist(),
                   "residual": max(disc.residuals), "clusters": errors,
                   "seconds": time.time() - start})
      print(f"hdg  p={degree} n={n:4d} unknowns={disc.n_unknowns():8d} " + "  ".join(
        f"{name}: {errors[name]['eigenvalue_error']:.3e} {errors[name]['gap']:.3e}"
        for name in CLUSTERS), flush=True)
  return runs


def thin_study(quick):
  exponents = range(3, 5) if quick else range(3, 11)
  reference = dict(THIN_REFERENCE, degree=3) if quick else THIN_REFERENCE
  runs = []
  for exponent in exponents:
    eps = 2. ** -exponent
    start = time.time()
    sectors = thin.thin_eigenpairs(eps, **reference)
    errors = thin.errors(eps, sectors, reference["degree"])
    values = sorted(float(v) for k in (0, 1) for v in sectors[k][2])
    run = {"eps": eps, "unknowns": {k: sectors[k][1].n_unknowns() for k in (0, 1)},
           "eigenvalues": values[:4], "order_as_limit": thin.global_order_ok(sectors),
           "residual": max(max(sectors[k][1].residuals) for k in (0, 1)), "clusters": errors}
    if not quick:
      check = thin.thin_eigenpairs(eps, **THIN_CHECK)
      check_values = sorted(float(v) for k in (0, 1) for v in check[k][2])[:4]
      run["check_eigenvalues"] = check_values
      run["check_difference"] = float(np.max(np.abs(np.array(check_values) - values[:4])))
    run["seconds"] = time.time() - start
    runs.append(run)
    print(f"thin eps=2^-{exponent:<2d} " + "  ".join(
      f"{name}: {errors[name]['eigenvalue_error']:.4e} {errors[name]['gap']:.4e}"
      for name in CLUSTERS) + (f"  check: {run['check_difference']:.1e}" if not quick else ""),
      flush=True)
  return runs


def combined_study(hdg_runs, thin_runs, degree=1):
  """max_{1 <= j <= 4} |lambda_{eps,j} - lambda_{h,j}| for all eps and h."""
  table = []
  for thin_run in thin_runs:
    row = []
    for run in [r for r in hdg_runs if r["degree"] == degree]:
      diff = max(abs(a - b) for a, b in zip(thin_run["eigenvalues"], run["eigenvalues"][:4]))
      row.append({"n": run["n"], "h": run["h"], "difference": diff})
    table.append({"eps": thin_run["eps"], "row": row})
  return table


# --------------------------------------------------------------------------------------------------
# LaTeX tables
# --------------------------------------------------------------------------------------------------

def sci(value):
  """1.23e-04 -> 1.23\\cdot10^{-4} (for math mode)."""
  mantissa, exponent = f"{value:.2e}".split("e")
  return mantissa if int(exponent) == 0 else f"{mantissa}\\cdot10^{{{int(exponent)}}}"


def order(value):
  return "--" if value is None else f"{value:.2f}"


def eps_label(eps):
  return f"$2^{{-{round(-math.log2(eps))}}}$"


def thin_tables(runs):
  eps = [run["eps"] for run in runs]
  lines = []
  for quantity, label in [("eigenvalue_error", "e"), ("gap", "g")]:
    columns = []
    for name in CLUSTERS:
      values = [run["clusters"][name][quantity] for run in runs]
      columns.append((values, common.eoc(values, eps)))
    body = []
    for i, run in enumerate(runs):
      cells = [eps_label(run["eps"])]
      for values, orders in columns:
        cells += [f"${sci(values[i])}$", order(orders[i])]
      body.append(" & ".join(cells) + r" \\")
    lines.append((quantity, body))
  return lines


def write_tables(hdg_runs, thin_runs, combined):
  os.makedirs(RESULTS, exist_ok=True)
  # thin domain
  if thin_runs:
    for quantity, body in thin_tables(thin_runs):
      with open(os.path.join(RESULTS, f"table_thin_{quantity}.tex"), "w") as file:
        file.write("\n".join(body) + "\n")
  # HDG
  if hdg_runs:
    for quantity in ["eigenvalue_error", "gap"]:
      body = []
      for degree in sorted(set(run["degree"] for run in hdg_runs)):
        runs = [run for run in hdg_runs if run["degree"] == degree]
        sizes = [run["h"] for run in runs]
        columns = []
        for name in CLUSTERS:
          values = [run["clusters"][name][quantity] for run in runs]
          columns.append((values, common.eoc(values, sizes)))
        for i, run in enumerate(runs):
          cells = [f"\\multirow{{{len(runs)}}}{{*}}{{{degree}}}" if i == 0 else "",
                   str(run["n"]), str(run["unknowns"])]
          for values, orders in columns:
            cells += [f"${sci(values[i])}$", order(orders[i])]
          body.append(" & ".join(cells) + r" \\")
        body.append(r"\midrule")
      body = body[:-1]
      with open(os.path.join(RESULTS, f"table_hdg_{quantity}.tex"), "w") as file:
        file.write("\n".join(body) + "\n")
  # combined
  if combined:
    body = []
    for row in combined:
      cells = [eps_label(row["eps"])] + [f"${sci(entry['difference'])}$" for entry in row["row"]]
      body.append(" & ".join(cells) + r" \\")
    with open(os.path.join(RESULTS, "table_combined.tex"), "w") as file:
      file.write("\n".join(body) + "\n")


def main():
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--quick", action="store_true", help="coarse levels only (a few seconds)")
  parser.add_argument("--only", choices=["hdg", "thin", "combined"], nargs="+",
                      default=["hdg", "thin", "combined"])
  parser.add_argument("--tables", action="store_true",
                      help="only write the tables from the results of a previous run")
  options = parser.parse_args()
  if options.tables:
    results = {}
    for name in ["hdg", "thin", "combined"]:
      with open(os.path.join(RESULTS, f"{name}.json")) as file:
        results[name] = json.load(file)["runs"]
    write_tables(results["hdg"], results["thin"], results["combined"])
    return

  print("Limit eigenvalues:", ", ".join(f"{name} = {limit.eigenvalue(cluster):.15f}" for name,
                                        cluster in limit.clusters().items()), flush=True)
  os.makedirs(RESULTS, exist_ok=True)
  suffix = "_quick" if options.quick else ""
  meta = {"python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
          "numpy": np.__version__, "scipy": __import__("scipy").__version__,
          "limit_eigenvalues": {name: limit.eigenvalue(cluster)
                                for name, cluster in limit.clusters().items()}}

  def path(name):
    return os.path.join(RESULTS, f"{name}{suffix}.json")

  def load(name):
    with open(path(name)) as file:
      return json.load(file)["runs"]

  hdg_runs = thin_runs = combined = None
  if "hdg" in options.only:
    hdg_runs = hdg_study(options.quick)
    with open(path("hdg"), "w") as file:
      json.dump({"meta": meta, "runs": hdg_runs}, file, indent=1)
  if "thin" in options.only:
    thin_runs = thin_study(options.quick)
    with open(path("thin"), "w") as file:
      json.dump({"meta": dict(meta, reference=THIN_REFERENCE, check=THIN_CHECK),
                 "runs": thin_runs}, file, indent=1)
  if "combined" in options.only:
    hdg_runs = hdg_runs if hdg_runs is not None else load("hdg")
    thin_runs = thin_runs if thin_runs is not None else load("thin")
    combined = combined_study(hdg_runs, thin_runs)
    with open(path("combined"), "w") as file:
      json.dump({"meta": meta, "runs": combined}, file, indent=1)
  if not options.quick:
    write_tables(hdg_runs, thin_runs, combined)


if __name__ == "__main__":
  main()
