#!/usr/bin/env python3

import sys
import json
import argparse
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="name of x variable")
parser.add_argument("-y", help="name of y variable(s)")
parser.add_argument("--group-by", help="how to group input data to display")
parser.add_argument("--xlabel", help="label of x-axis, default=name of x variable")
parser.add_argument("--ylabel", help="label of y-axis, default=name of y variable")
parser.add_argument("--save", help="save plot instead of showing it")
parser.add_argument("--title", help="title of plot")

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

args = parser.parse_args()

objs = []

for line in sys.stdin:
    try:
        obj = json.loads(line)
        objs.append(obj)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e} in line: {line.strip()}", file=sys.stderr)

objs = sorted(objs, key=lambda o: o[args.x])

if args.group_by:
    group_variable = args.group_by
    groups = defaultdict(list)
    for o in objs:
        groups[o[group_variable]].append(o)
else:
    groups = {"group": objs}

for group_name, group in groups.items():
    xys = np.array([(o[args.x],o[args.y]) for o in group])
    plt.plot(xys[:, 0], xys[:, 1], label=group_name)

plt.xlabel(args.xlabel or args.x)
plt.ylabel(args.ylabel or args.y[0])
plt.title(args.title)
plt.legend()

if args.save:
    plt.savefig(args.save)
else:
    plt.show()
