import numpy as np
import sys
import time
import argparse
import pandas
import h5py
import yaml
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)


parser = argparse.ArgumentParser(description="histogram by Joseph Holten")
parser.add_argument("-i", "--input", help="input domain file", required=True)
parser.add_argument("-o", "--output", help="output file", default="histogram.png")
parser.add_argument("-b", "--bins", help="number of histogram bins", type=int, default=10)
parser.add_argument("--title", default="")
parser.add_argument("--cdf", action="store_true", help="plot cumulative density instead of mass density")
parser.add_argument("--csv", action="store_true", help="emit histogram as CSV on stdout (pipe into plot.py) instead of plotting")
args = parser.parse_args()

with h5py.File(args.input, "r") as f:
    points = f["domain/points"][:]
    edges  = f["domain/edges"][:]
    unit_length = f["domain"].attrs.get("unit_length", "")
    props = f["domain/properties"][:]

endpoints = points[edges]

he = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)
virtual = props[:, 0] == 0.
moment2 = w1*w2**3/12
Cq = props[:, 1:7]
cz = props[:, 0]
nz = ~virtual                      # cz != 0
ce = np.full(len(props), np.nan)   # or np.nan, your choice for virtual edges
ce = np.max(np.sqrt(Cq), axis=1)
ts = he / ce                       # virtual -> 0 if ce=inf, nan if ce=nan

bins = np.logspace(np.log10(ts.min()), np.log10(ts.max()), args.bins)

if args.csv:
    counts, edges = np.histogram(ts, bins=bins, weights=np.ones(len(ts))/len(ts))
    if args.cdf: counts = np.cumsum(counts)
    centers = np.sqrt(edges[:-1] * edges[1:])  # geometric centers of log bins
    mode = centers[np.argmax(counts)]
    col = "cdf" if args.cdf else "density"
    sys.stdout.write(f"# min={ts.min():.6g} median={np.median(ts):.6g} mode={mode:.6g} max={ts.max():.6g}\n")
    pandas.DataFrame({"ts": centers, col: counts}).to_csv(sys.stdout, index=False)
    sys.exit(0)

counts, edges, _ = plt.hist(ts, bins=bins,
                            weights=np.ones(len(ts))/len(ts), cumulative=args.cdf)
plt.xscale("log")
if not args.cdf: plt.yscale("log")

mode = np.sqrt(edges[np.argmax(counts)] * edges[np.argmax(counts) + 1])  # geometric center of tallest bin
stats = {
    "min":    (ts.min(),       "C1"),
    "median": (np.median(ts),  "C2"),
    "mode":   (mode,           "C3"),
    "max":    (ts.max(),       "C4"),
}
for label, (val, color) in stats.items():
    plt.axvline(val, color=color, linestyle="--", label=f"{label} = {val:.3g}")
plt.legend()

plt.xlabel("timestep ts")
plt.ylabel("cumulative density" if args.cdf else "density")
if args.title: plt.title(args.title)
if args.output: plt.savefig(args.output)
plt.show()
