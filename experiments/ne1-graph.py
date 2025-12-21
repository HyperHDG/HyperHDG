#!/usr/bin/env python3

import sys
import argparse
import matplotlib.pyplot as plt
import pandas as pd

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
parser.add_argument("-f", "--format", help="format of input, csv|json")
parser.add_argument("--lines", help="new line delimited json", default=True, action="store_true")

args = parser.parse_args()

plot_func = plt.plot if not args.scatter else plt.scatter

match args.format:
    case "csv": df = pd.read_csv(sys.stdin, comment="#")
    case "json": df = pd.read_json(sys.stdin, lines=args.lines)
    case _:
        print("ERROR: unrecognized format", args.format, file=sys.stderr)
        sys.exit(1)

for name, group in df.groupby(args.group_by) if args.group_by else [(None,df)]:
    plot_func(group[args.x], group[args.y], label=f"{args.group_by}={name}", marker="+")

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
