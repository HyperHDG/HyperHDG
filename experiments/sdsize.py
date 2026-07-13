#!/usr/bin/env python3
"""
sdsize.py — distribution of net2as subdomain sizes vs the per-axis count p.

Reads newline-delimited JSON from stdin, one record per p:

    {"p": <int>, "sizes": [<nonzero subdomain DoF counts>], "nidle": <int>}

and draws a box plot of the (nonzero) subdomain sizes for each p on a log-y
axis. "nidle" is the number of size-0 entries net2as printed for that p; these
are NOT empty network patches but idle-rank placeholders: net2as_alloc_ds keeps
max(1, sz) slots per rank, so when there are fewer subdomains than MPI ranks
every idle rank still emits one size-0 line. They are annotated above each box
so the count is visible but they don't distort the distribution.

Example
-------
    jq ... log.json | sdsize.py --save sizes.png

Author
------
Joseph Holten, KIT, 2026.
"""

import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="net2as subdomain size distribution")
parser.add_argument("--save", help="save plot instead of showing it")
parser.add_argument("--title", default="subdomain size distribution")
parser.add_argument("--figsize", default="6,6", help="figure size 'w,h' in inches")
args = parser.parse_args()

recs = [json.loads(line) for line in sys.stdin if line.strip()]
recs.sort(key=lambda r: r["p"])

ps = [r["p"] for r in recs]
sizes = [np.asarray(r["sizes"], dtype=float) for r in recs]
nidle = [r.get("nidle", 0) for r in recs]
pos = np.arange(len(ps))  # categorical x: even spacing for geometric p

w, h = (float(s) for s in args.figsize.split(","))
plt.figure(figsize=(w, h))

plt.boxplot(sizes, positions=pos, widths=0.6, whis=(0, 100),
            showmeans=True, meanline=True,
            medianprops=dict(color="C0"), meanprops=dict(color="C1", ls="--"))

# annotate the idle-rank placeholder count above each box (only where nonzero)
for x, s, ni in zip(pos, sizes, nidle):
    if ni:
        plt.text(x, s.max() * 1.6, f"+{ni}\nidle", ha="center", va="bottom",
                 fontsize=8, color="0.4")

plt.yscale("log")
plt.xticks(pos, ps)
plt.xlabel("subdomains per axis  $p$  (smaller $H$ →)")
plt.ylabel("subdomain size  [DoFs]")
plt.title(args.title)
plt.margins(y=0.15)
plt.grid(axis="y", which="both", ls=":", alpha=0.5)

# legend for the median/mean lines
from matplotlib.lines import Line2D
plt.legend(handles=[Line2D([], [], color="C0", label="median"),
                    Line2D([], [], color="C1", ls="--", label="mean"),
                    Line2D([], [], color="0.4", label="whiskers: min/max")],
           loc="best")

plt.tight_layout()
if args.save:
    plt.savefig(args.save, bbox_inches="tight", pad_inches=.05)
    print(f"saved {args.save}")
else:
    plt.show()
