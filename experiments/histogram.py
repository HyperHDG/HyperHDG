import numpy as np
import time
import argparse
import pandas
import h5py
import yaml
import matplotlib.pyplot as plt

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)


parser = argparse.ArgumentParser(description="histogram by Joseph Holten")
parser.add_argument("-i", "--input", help="input domain file", required=True)
parser.add_argument("-o", "--output", help="output file", default="histogram.png")
parser.add_argument("-b", "--bins", help="number of histogram bins", type=int, default=10)
parser.add_argument("--title", default="")
args = parser.parse_args()

with h5py.File(args.input, "r") as f:
    points = f["domain/points"][:]
    edges  = f["domain/edges"][:]
    unit_length = f["domain"].attrs["unit_length"]

endpoints = points[edges]
lengths = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)

bins = np.logspace(np.log10(lengths.min()), np.log10(lengths.max()), args.bins + 1)
plt.hist(lengths, bins=bins, weights=np.ones_like(lengths) / len(lengths))
plt.xscale("log")
plt.yscale("log")
plt.xlabel(f"edge length {unit_length}")
plt.ylabel("density")
plt.title(args.title)
if args.output: plt.savefig(args.output)
plt.show()
