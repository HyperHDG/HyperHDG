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
parser.add_argument("-b", "--bins", help="number of histogram bins", type=int, default=100)
parser.add_argument("--title", default="")
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
print("fraction massless edges: ", sum(virtual)/len(he))

bins = np.logspace(np.log10(he.min()), np.log10(he.max()), args.bins)
nv, v = he[~virtual], he[virtual]
N = len(he)  # total edges, so the two classes are normalized to the same denominator
plt.hist(nv, bins=bins, label=r"$m_e>0$", weights=np.ones(len(nv))/N,
         histtype="step", linewidth=2, cumulative=True)
plt.hist(v,  bins=bins, label=r"$m_e=0$", weights=np.ones(len(v))/N,
         histtype="step", linewidth=2, cumulative=True)
plt.xscale("log")
plt.yscale("log")
plt.title("cumulative density of edge lengths")
plt.xlabel("$h_e$")
plt.ylabel("density")
plt.legend()
plt.show()

props = props[~virtual, :]
he = he[~virtual]

mass = props[:, 0]
w1 = props[:, 13]
w2 = props[:, 14]

moment1 = w1**3*w2/12
moment2 = w1*w2**3/12

density = mass/(w1*w2*he)

Cq = props[:, 1:7]
Cz = np.stack((mass/he, mass/he, mass/he,
                density * (moment1+moment2), density*moment1, density*moment2), axis=1)
ce = np.max(np.sqrt(Cq), axis=1)
ts = he / ce

bins = np.logspace(np.log10(ts.min()), np.log10(ts.max()), args.bins)

if args.csv:
    counts, edges = np.histogram(ts, bins=bins, weights=np.ones(len(ts))/len(ts))
    centers = np.sqrt(edges[:-1] * edges[1:])  # geometric centers of log bins
    mode = centers[np.argmax(counts)]
    sys.stdout.write(f"# min={ts.min():.6g} median={np.median(ts):.6g} mode={mode:.6g} max={ts.max():.6g}\n")
    pandas.DataFrame({"ts": centers, col: counts}).to_csv(sys.stdout, index=False)
    sys.exit(0)

counts, edges, _ = plt.hist(ts, bins=bins,
                            weights=np.ones(len(ts))/len(ts))
plt.xscale("log")
plt.yscale("log")

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
plt.ylabel("density")
if args.title: plt.title(args.title)
if args.output: plt.savefig(args.output)
plt.show()
