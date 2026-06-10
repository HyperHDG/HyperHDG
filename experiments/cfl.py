import numpy as np
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
args = parser.parse_args()

with h5py.File(args.input, "r") as f:
    points = f["domain/points"][:]
    edges  = f["domain/edges"][:]
    unit_length = f["domain"].attrs.get("unit_length", "")
    props = f["domain/properties"][:]

endpoints = points[edges]

he = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)
virtual = props[:, 0] == 0.
Cq = props[:, 1:7]
cz = props[:, 0]
nz = ~virtual                      # cz != 0
ce = np.full(len(props), np.nan)   # or np.nan, your choice for virtual edges
ce = np.max(np.sqrt(Cq), axis=1)
ts = he / ce                       # virtual -> 0 if ce=inf, nan if ce=nan

# TODO
counts, edges, _ = plt.hist(ts, bins=np.logspace(np.log10(ts.min()), np.log10(ts.max()), args.bins),
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
