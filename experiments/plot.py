#!/usr/bin/env python3

import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import numpy as np


def fmt_names(names):
  if names is None: return ''
  return ",".join(map(str, names))


def plt_legend2(legend_title=None, lbbox=None):
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.insert(0, Patch(color="none", visible=False))
    labels.insert(0, legend_title)

    lbbox = lbbox.split(';')
    bbox = tuple(map(float,lbbox[1].split(','))) if len(lbbox) > 1 else None

    plt.legend(handles, labels, bbox_to_anchor=bbox, loc=lbbox[0])


def reference_triangle_loglog(rate, x0, y0, tx, ty, **kw):
    x = tx(np.array(x0))
    y = x**rate
    y /= y[0]
    y *= y0
    y = ty(y)

    # up
    xs = [x[0], x[1], x[1], x[0]]
    ys = [y[0], y[0], y[1], y[0]]

    # TODO: if triangle is upside down, make the text be on the left instead of the right side
    plt.plot(xs, ys, **kw)
    plt.text(np.sqrt(x[0]*x[1]), ys[0]*.9, '1', ha='center', va='top')
    plt.text(xs[1]*1.05, np.sqrt(y[0]*y[1]), str(rate), ha='left', va='center')

    # down
    # swap the 0th and 1st elements of x and y, use max(*y) and min(*x) instead and va=bottom, ha=right instead

#### MAIN ####

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="name of x variable")
parser.add_argument("-y", help="name of y variable(s)")
parser.add_argument("-g", "--group-by", help="how to group input data to display")
parser.add_argument("--xlabel", help="label of x-axis, default=name of x variable")
parser.add_argument("--ylabel", help="label of y-axis, default=name of y variable")
parser.add_argument("--xbase", help="base of x-axis", type=int, default=10)
parser.add_argument("--ybase", help="base of y-axis", type=int, default=10)
parser.add_argument("--save", help="save plot instead of showing it")
parser.add_argument("--title", help="title of plot")
parser.add_argument("--log", help="axis to apply log scale")
parser.add_argument("--scatter", help="show as scatter plot", action="store_true")
parser.add_argument("--nshow", help="don't show the plot", action="store_true")
parser.add_argument("-f", "--format", help="format of input, csv|json")
parser.add_argument("--lines", help="new line delimited json", default=True, action="store_true")
parser.add_argument("-w", "--where", help="select subset")
parser.add_argument("--legend", help="legend loc '<loc: str>;<bbox: float,float>'", default='best')
parser.add_argument("--trans", help="transform input 'f(x),g(y)'")
parser.add_argument("--ref", help="generate reference triangle 'rate;x0,x1;y0'")

args = parser.parse_args()

plot_func = plt.plot if not args.scatter else plt.scatter
match args.format:
    case "csv": df = pd.read_csv(sys.stdin, comment="#")
    case "json": df = pd.read_json(sys.stdin, lines=args.lines)
    case _:
        print("ERROR: unrecognized format", args.format, file=sys.stderr)
        sys.exit(1)

if args.where:
    df = df.query(args.where)

if args.trans:
    ts = args.trans.split(',')
    assert len(ts) == 2, "--trans format 'f(x),g(y)'"
    tx, ty = [eval(f"lambda {v}: {t}") for t, v in zip(ts, 'xy')]
else:
    tx, ty = lambda x: x, lambda y: y

for names, group in df.groupby(args.group_by.split(',')) if args.group_by else [(None,df)]:
    plot_func(tx(group[args.x]), ty(group[args.y]), label=fmt_names(names), marker="+")

if args.ref:
    rate, x, y0 = args.ref.split(';')
    reference_triangle_loglog(int(rate), [float(xi) for xi in x.split(',')],
      float(y0), tx, ty, linewidth=1, color='.5')

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
if args.log:
    if "x" in args.log: plt.xscale("log", base=args.xbase)
    if "y" in args.log: plt.yscale("log", base=args.ybase)
plt.xlabel(args.xlabel or args.x)
plt.ylabel(args.ylabel or args.y)
plt.title(args.title)
if args.group_by and args.legend: plt_legend2(legend_title=args.group_by, lbbox=args.legend)
plt.gca().set_box_aspect(1)
plt.tight_layout()

if args.save: plt.savefig(args.save, bbox_inches="tight", pad_inches=0)
if not args.nshow: plt.show()
