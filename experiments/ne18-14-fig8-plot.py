#!/usr/bin/env python
"""Plot the gortz.pdf figure-8 reproduction: per-iteration K-norm error of PCG with the
net2as subspace decomposition preconditioner, one curve per coarse scale H = 1/p, from
the newline-delimited JSON log written by ne18-14-fig8-grid.sh (network -err_monitor).
Also prints Table-2-style (average, worst) convergence rates tau_(l) = e_l / e_{l-1}
for l >= 2 (gortz.pdf eq. 6.2)."""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("log", help="newline-delimited JSON, one run per line")
    ap.add_argument("-o", "--output", default="fig8.png")
    args = ap.parse_args()

    with open(args.log) as f:
        runs = [json.loads(line) for line in f if line.strip()]

    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    styles = ["-", "--", "-.", ":"]
    print(f"{'H':>8} {'iters':>6} {'cond':>10} {'tau_avg':>8} {'tau_max':>8}")
    for k, run in enumerate(runs):
        p = run.get("Hinv", run.get("net2as_p", k))
        mon = run["ksp_err_monitor"]
        it = np.array([m["it"] for m in mon])
        err = np.array([m["err_K"] for m in mon])
        # rates tau_(l) for l >= 2, excluding iterates saturated at the
        # reference-solution accuracy
        valid = err > err[0] * 1e-12
        rates = err[1:][valid[1:]] / err[:-1][valid[1:]]
        rates = rates[1:]
        print(f"{'1/' + str(p):>8} {it[-1]:>6} {run.get('cond', float('nan')):>10.4g} "
              f"{rates.mean():>8.3f} {rates.max():>8.3f}")
        ax.semilogy(it, err, styles[k % len(styles)], label=f"$H = 1/{p}$")

    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$\|u - u^{(\ell)}\|_K$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print("saved", args.output)


if __name__ == "__main__":
    main()
