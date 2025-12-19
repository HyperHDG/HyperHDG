#!/usr/bin/env python3

import sys
import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="name of x variable")
parser.add_argument("-y", help="name of y variable(s)")
parser.add_argument("-g", "--group-by", help="how to group input data to display")
parser.add_argument("--xlabel", help="label of x-axis, default=name of x variable")
parser.add_argument("--ylabel", help="label of y-axis, default=name of y variable")
parser.add_argument("--save", help="save plot instead of showing it")
parser.add_argument("--title", help="title of plot")
parser.add_argument("--log", help="axis to apply log scale")
parser.add_argument("--scatter", help="show as scatter plot", action="store_true")
parser.add_argument("--nshow", help="don't show the plot", action="store_true")

args = parser.parse_args()

plot_func = plt.plot if not args.scatter else plt.scatter

objs = (json.loads(line) for line in sys.stdin)

if args.group_by:
    group_variable = args.group_by
    groups = defaultdict(list)
    for o in objs:
        groups[o[group_variable]].append(o)
else:
    groups = {"default_group": objs}

sorted_groups = sorted(groups.items(), key=lambda kv: kv[0])

for group_name, group in sorted_groups:
    xys = np.array([(o[args.x],o[args.y]) for o in group])
    plot_func(xys[:, 0], xys[:, 1], label=f"{args.group_by}={group_name}", marker="+")

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
if args.log:
    if "x" in args.log:
        plt.xscale("log")
    if "y" in args.log:
        plt.yscale("log")
plt.xlabel(args.xlabel or args.x)
plt.ylabel(args.ylabel or args.y)
plt.title(args.title)
if args.group_by:
    plt.legend()

if args.save: plt.savefig(args.save)
if not args.nshow: plt.show()
