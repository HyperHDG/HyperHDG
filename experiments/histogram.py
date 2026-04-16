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
args = parser.parse_args()

with h5py.File(args.input, "r") as f:
    points = f["domain/points"][:]
    edges  = f["domain/edges"][:]
    unit_length = f["domain"].attrs["unit_length"]
    props = f["domain/properties"][:]

endpoints = points[edges]
lengths = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)

kG1A, kG2A = props[:, 1], props[:, 2]
E1I1, E2I2 = props[:, 4], props[:, 5]
k_shear = np.sqrt(np.maximum(kG1A / E1I1, kG2A / E2I2))
ratio = k_shear * lengths

lo = min(lengths.min(), (1/k_shear).min())
hi = max(lengths.max(), (1/k_shear).max())

plt.hist2d(np.log10(lengths), np.log10(k_shear),
           bins=args.bins, weights=np.ones(len(lengths))/len(lengths),
           cmin=1e-7, norm=LogNorm())
plt.colorbar(label="density")
plt.xlabel(f"log10(L) {unit_length}")
plt.ylabel("log10($k_{shear}$)")
if args.output: plt.savefig(args.output)
plt.show()
