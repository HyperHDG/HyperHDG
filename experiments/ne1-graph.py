#!/usr/bin/env python3

import sys
import json
import argparse
from pprint import pprint
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="name of x variable")
parser.add_argument("-y", help="name of y variable")
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

plt.xlabel(args.x)
plt.xlabel(args.y)

xs = [o[args.x] for o in objs]
ys = [o[args.y] for o in objs]

plt.plot(xs, ys)

if args.save:
    plt.savefig(args.save)
else:
    plt.show()
