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
    props = f["domain/properties"][:]

endpoints = points[edges]
lengths = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)

kG1A, kG2A = props[:, 1], props[:, 2]
E1I1, E2I2 = props[:, 4], props[:, 5]
k_shear = np.sqrt(np.maximum(kG1A / E1I1, kG2A / E2I2))

lo = min(lengths.min(), (1/k_shear).min())
hi = max(lengths.max(), (1/k_shear).max())
bins = np.logspace(np.log10(lengths.min()), np.log10(lengths.max()), args.bins + 1)

plt.hist(lengths, bins=bins, weights=np.ones_like(lengths) / len(lengths), histtype="step", label="lengths")
plt.hist(1/k_shear, bins=bins, weights=np.ones_like(k_shear)/len(k_shear),
         histtype="step", label=r"$k_shear=\sqrt{EI/kGA}$")

kG1A, kG2A = props[:, 1], props[:, 2]
E1I1, E2I2 = props[:, 4], props[:, 5]
k_shear = np.sqrt(np.maximum(kG1A / E1I1, kG2A / E2I2))

plt.xscale("log")
plt.yscale("log")
plt.xlabel(f"edge length {unit_length}")
plt.ylabel("density")
plt.title(args.title)
if args.output: plt.savefig(args.output)
plt.show()
