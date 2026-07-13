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
parser.add_argument("--cdf", help="make timestep histogram cdf", action="store_true")
parser.add_argument("--noshow", action="store_true")
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
print("# fraction massless edges: ", sum(virtual)/len(he))

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
if not args.noshow: plt.show()

# 2d spatial distribution of areal mass density: each edge's mass sits at its
# midpoint, summed per xy-bin and divided by the bin area -> mass per area
mid = endpoints.mean(axis=1)            # (N, space_dim) edge midpoints
mass_per_bin, xe, ye = np.histogram2d(mid[:, 0], mid[:, 1], bins=args.bins, weights=props[:, 0])
bin_area = np.diff(xe)[:, None] * np.diff(ye)[None, :]
areal_density = mass_per_bin / bin_area
rel = areal_density / areal_density.mean() - 1  # fractional deviation from mean (gsm factor cancels)
lim = np.abs(rel).max()
# .T: histogram2d is [x, y]-indexed, pcolormesh wants [row=y, col=x]
plt.pcolormesh(xe, ye, 100 * rel.T, cmap="RdBu_r", vmin=-100 * lim, vmax=100 * lim)
plt.colorbar(label="areal weight rel. to mean [%]")
plt.gca().set_aspect("equal")
plt.title(f"areal mass density, mean = {areal_density.mean()*1e15:.2f} gsm")
plt.xlabel(f"x {unit_length}")
plt.ylabel(f"y {unit_length}")
if not args.noshow: plt.show()

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
ce = np.max(np.sqrt(Cq/Cz), axis=1)
ts = he / ce

bins = np.logspace(np.log10(ce.min()), np.log10(ce.max()), args.bins)
plt.hist(ce, bins=bins, label=r"$c_e$", weights=np.ones(len(ce))/len(ce),
         histtype="step", linewidth=2)
plt.title("density of wavespeed $c_e$ across varying edges $e$")
plt.ylabel("density")
plt.xlabel("$c_e$")
plt.yscale("log")
if not args.noshow: plt.show()

bins = np.logspace(np.log10(ts.min()), np.log10(ts.max()), args.bins)

if args.csv:
    counts, edges = np.histogram(ts, bins=bins, weights=np.ones(len(ts))/len(ts))
    if args.cdf:
        counts = np.cumsum(counts)  # np.histogram has no cumulative flag; prefix-sum the counts
    centers = np.sqrt(edges[:-1] * edges[1:])  # geometric centers of log bins
    mode = centers[np.argmax(counts)]
    col = "cdf" if args.cdf else "density"
    sys.stdout.write(f"# min={ts.min():.6g} median={np.median(ts):.6g} mode={mode:.6g} max={ts.max():.6g}\n")
    pandas.DataFrame({"ts": centers, col: counts}).to_csv(sys.stdout, index=False)
    sys.exit(0)

counts, edges, _ = plt.hist(ts, bins=bins,
                            weights=np.ones(len(ts))/len(ts), cumulative=args.cdf)
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
plt.ylabel(("cumulative " if args.cdf else "") + "density")
if args.title: plt.title(args.title)
if args.output: plt.savefig(args.output)
if not args.noshow: plt.show()
