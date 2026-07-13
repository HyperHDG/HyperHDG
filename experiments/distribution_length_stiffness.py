import numpy as np
import time
import argparse
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def tprint(*args, **kwargs):
    print(f"[{time.strftime('%H:%M:%S')}]", *args, **kwargs)

parser = argparse.ArgumentParser(description="length-vs-stiffness histogram by Joseph Holten")
parser.add_argument("-i", "--input", help="input domain file", required=True)
parser.add_argument("-o", "--output", help="output file", default="histogram.png")
parser.add_argument("-b", "--bins", help="number of histogram bins", type=int, default=10)
parser.add_argument("-c", "--component", help=f"stiffness component index",
                    type=int, default=0)
parser.add_argument("--title", default="")
parser.add_argument("--no-show", action="store_true")
args = parser.parse_args()

with h5py.File(args.input, "r") as f:
    points = f["domain/points"][:]
    edges  = f["domain/edges"][:]
    unit_length = f["domain"].attrs.get("unit_length", "")
    props = f["domain/properties"][:]

endpoints = points[edges]
lengths = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1)

stiffness = props[:, args.component]

plt.hist2d(np.log10(lengths), np.log10(stiffness),
           bins=args.bins, weights=np.ones(len(lengths))/len(lengths),
           cmin=1e-7, norm=LogNorm())
plt.colorbar(label="density")
plt.xlabel(f"log10(L) {unit_length}")
plt.ylabel(f"log10(props_{args.component})")
if args.title: plt.title(args.title)
if args.output: plt.savefig(args.output)
if not args.no_show: plt.show()
