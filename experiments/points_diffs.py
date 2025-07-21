#!/usr/bin/env python3

import numpy as np
import sys

if len(sys.argv) < 2:
    print(f"ERROR: usage: {sys.argv[0]} <a.txt> <b.txt>")

a = np.loadtxt(sys.argv[1])
b = np.loadtxt(sys.argv[2])

diffs = np.linalg.norm(a-b, axis=1)
print(f"avg diff={np.mean(diffs)}, max diff={np.max(diffs)}")
