#!/usr/bin/env python3

import sys
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path


def fmt_names(names):
    def fmt(x):
        match x:
            case float() | np.floating(): return f"{x:.3e}"
            case _: return str(x)
    return ",".join(map(fmt, names))


def tex_escape(s):
    # escape special chars in column-name-derived strings (user-provided labels are passed through untouched)
    return str(s).replace('\\', '\\textbackslash{}').replace('_', r'\_').replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')


LEGEND_POS_MAP = {
    'best': 'outer north east',
    'upper right': 'north east',
    'upper left': 'north west',
    'lower right': 'south east',
    'lower left': 'south west',
    'center right': 'east',
    'center left': 'west',
    'upper center': 'north',
    'lower center': 'south',
    'center': 'center',
}


def collect_plots(df, args, tx, ty):
    """Iterate the data once, yielding (idx, name0, series, refs) per outer-group.

    series: list of (xs, ys, label) — already sorted, transformed, EOC-computed
    refs:   list of (rate, xt, yt) — reference triangle vertices in plot coords
    """
    refs = []
    if args.ref and not args.eoc:
        for ref in args.ref.split('|'):
            rate_s, xs_s, y0_s = ref.split(';')
            rate = float(rate_s)
            xs_in = np.array([float(xi) for xi in xs_s.split(',')])
            y0 = float(y0_s)
            xt = tx(xs_in)
            yt = xt**rate
            yt = yt / yt[0] * y0
            yt = ty(yt)
            refs.append((rate, xt, yt))

    outer = enumerate(df.groupby(args.group0)) if args.group0 else [(0, (None, df))]
    for idx, (name0, df0) in outer:
        series = []
        inner = df0.groupby(args.group_by.split(',')) if args.group_by else [("", df0)]
        for names, group in inner:
            if name0 is not None: names = [name0] + list(names)
            sgroup = group[[args.x, args.y]].sort_values(args.x)
            xs, ys = tx(sgroup[args.x].to_numpy()), ty(sgroup[args.y].to_numpy())
            if args.eoc:
                xs, ys = xs[1:], np.log(ys[1:]/ys[:-1]) / np.log(xs[1:]/xs[:-1])
            label = fmt_names(names) if (args.group_by or name0 is not None) else ""
            series.append((xs, ys, label))
        yield idx, name0, series, refs


def save_path(args, idx, name0, base):
    """Apply group0 suffix to a save path if needed."""
    p = Path(base)
    if name0 is None: return p
    return (p.parent / f"{p.stem}_{args.group0}{idx}").with_suffix(p.suffix)


def legend_title(args):
    return f"{args.group0},{args.group_by}" if args.group0 is not None else args.group_by


def draw_ref_mpl(rate, xt, yt, **kw):
    rate_lbl = int(rate) if rate.is_integer() else rate
    xs = [xt[0], xt[1], xt[1], xt[0]]
    ys = [yt[0], yt[0], yt[1], yt[0]]
    plt.plot(xs, ys, **kw)
    plt.text(np.sqrt(xt[0]*xt[1]), ys[0]*.9, '1', ha='center', va='top')
    plt.text(xs[1]*1.05, np.sqrt(yt[0]*yt[1]), str(rate_lbl), ha='left', va='center')


def render_matplotlib(idx, name0, series, refs, args, plot_func):
    w, h = map(float, args.figsize.split(','))
    plt.figure(figsize=(w, h))
    for xs, ys, label in series:
        plot_func(xs, ys, label=label, marker=args.marker)
    for rate, xt, yt in refs:
        draw_ref_mpl(rate, xt, yt, linewidth=1, color='.5')
    if args.comment:
        plt.text(1.0, -0.1, args.comment, transform=plt.gca().transAxes,
                 ha='right', va='top', fontsize=9, color=".5")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if args.log:
        if "x" in args.log: plt.xscale("log", base=args.xbase)
        if "y" in args.log: plt.yscale("log", base=args.ybase)
    plt.xlabel(args.xlabel or args.x)
    plt.ylabel(args.ylabel or args.y)
    plt.title(args.title)
    if args.group_by and args.legend:
        lbbox = args.legend.split(';')
        bbox = tuple(map(float, lbbox[1].split(','))) if len(lbbox) > 1 else None
        plt.legend(bbox_to_anchor=bbox, loc=lbbox[0], title=legend_title(args), alignment="left")
    plt.gca().set_box_aspect(1)
    plt.tight_layout()
    if args.save:
        for p in args.save.split(","):
            plt.savefig(save_path(args, idx, name0, p), bbox_inches="tight", pad_inches=.05)


def render_tikz(idx, name0, series, refs, args):
    prefix = Path(args.tikz)
    if str(prefix.parent) not in (".", ""):
        prefix.parent.mkdir(parents=True, exist_ok=True)

    tex_path = save_path(args, idx, name0, str(prefix) + ".tex")
    csv_path = save_path(args, idx, name0, str(prefix) + ".csv")

    frames = []
    for s_idx, (xs, ys, label) in enumerate(series):
        f = pd.DataFrame({args.x: xs, args.y: ys})
        f["series"] = s_idx
        f["legend entry"] = label
        frames.append(f)
    pd.concat(frames, ignore_index=True).to_csv(csv_path, index=False)

    legend_loc = args.legend.split(';')[0] if args.legend else 'best'
    legend_pos = LEGEND_POS_MAP.get(legend_loc, 'outer north east')

    L = [
        "% generated by experiments/plot.py",
        f"% standalone: latexmk -pdf {tex_path.name} (keep {csv_path.name} next to it)",
        "% manuscript: \\usepackage{standalone} \\usepackage{pgfplots} \\pgfplotsset{compat=newest},",
        "%   then \\input this file -- the preamble below is skipped, fonts come from the manuscript.",
        "%   \\renewcommand{\\figurewidth}{...} (default \\linewidth) and \\plotdatadir (default empty,",
        "%   set with trailing /) before \\input as needed.",
        "%   add \\usetikzlibrary{external} \\tikzexternalize to cache figures across manuscript compiles.",
        "\\documentclass{standalone}",
        "\\usepackage{amsmath}",
        "\\usepackage{pgfplots}",
        "\\pgfplotsset{compat=newest}",
        "\\begin{document}",
        "\\providecommand{\\figurewidth}{\\linewidth}",
        "\\providecommand{\\plotdatadir}{}",
        "\\begin{tikzpicture}",
    ]

    axis_opts = ["width=\\figurewidth"]
    if args.log:
        if "x" in args.log: axis_opts.append(f"xmode=log, log basis x={{{args.xbase}}}")
        if "y" in args.log: axis_opts.append(f"ymode=log, log basis y={{{args.ybase}}}")
    axis_opts.append(f"xlabel={{{args.xlabel or tex_escape(args.x)}}}")
    axis_opts.append(f"ylabel={{{args.ylabel or tex_escape(args.y)}}}")
    if args.title: axis_opts.append(f"title={{{args.title}}}")
    if args.group_by and args.legend:
        axis_opts.append(f"legend pos={legend_pos}")
        axis_opts.append(f"legend style={{title={{{tex_escape(legend_title(args))}}}, legend cell align=left}}")
        entries = ",\n    ".join("{" + tex_escape(lbl) + "}" for _, _, lbl in series)
        axis_opts.append(f"legend entries={{\n    {entries}\n  }}")

    L.extend([
        "\\begin{axis}[",
        ",\n".join("  " + o for o in axis_opts),
        "]",
    ])

    csv_name = csv_path.name
    L.extend([
        f"\\foreach \\i in {{0,...,{len(series)-1}}}{{",
        f"  \\addplot+[",
        f"    mark={args.marker},",
        f"    unbounded coords=discard,",
        f"    x filter/.expression={{\\thisrow{{series}} == \\i ? \\pgfmathresult : nan}},",
        f"  ] table[x={args.x}, y={args.y}, col sep=comma] {{\\plotdatadir {csv_name}}};",
        "}",
    ])

    for rate, xt, yt in refs:
        rate_lbl = int(rate) if rate.is_integer() else rate
        L.append(f"\\addplot[gray, thin, mark=none, forget plot] coordinates "
                 f"{{({xt[0]},{yt[0]}) ({xt[1]},{yt[0]}) ({xt[1]},{yt[1]}) ({xt[0]},{yt[0]})}};")
        xmid = (xt[0]*xt[1])**0.5
        ymid = (yt[0]*yt[1])**0.5
        L.append(f"\\node[below, gray, font=\\footnotesize] at (axis cs:{xmid},{yt[0]}) {{1}};")
        L.append(f"\\node[right, gray, font=\\footnotesize] at (axis cs:{xt[1]},{ymid}) {{{rate_lbl}}};")

    if args.comment:
        L.append(f"\\node[below right, font=\\scriptsize, gray] at (rel axis cs:1,0) {{{args.comment}}};")

    L.extend(["\\end{axis}", "\\end{tikzpicture}", "\\end{document}"])

    tex_path.write_text("\n".join(L) + "\n")
    print(f"wrote {tex_path} + {csv_path} ({len(series)} series)", file=sys.stderr)


#### MAIN ####

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="name of x variable")
parser.add_argument("-y", help="name of y variable(s)")
parser.add_argument("-g", "--group-by", help="how to group input data to display by color/ line")
parser.add_argument("--xlabel", help="label of x-axis, default=name of x variable")
parser.add_argument("--ylabel", help="label of y-axis, default=name of y variable")
parser.add_argument("--xbase", help="base of x-axis", type=int, default=10)
parser.add_argument("--ybase", help="base of y-axis", type=int, default=10)
parser.add_argument("--save", help="save plot instead of showing it")
parser.add_argument("--title", help="title of plot")
parser.add_argument("--log", help="axis to apply log scale")
parser.add_argument("--scatter", help="show as scatter plot", action="store_true")
parser.add_argument("--nshow", help="don't show the plot", action="store_true")
parser.add_argument("-f", "--format", help="format of input, csv|json", default="json")
parser.add_argument("--lines", help="new line delimited json", default=True, action="store_true")
parser.add_argument("-w", "--where", help="select subset")
parser.add_argument("--legend", help="legend loc '<loc: str>;<bbox: float,float>'", default='best')
parser.add_argument("--trans", help="transform input 'f(x),g(y)'")
parser.add_argument("--ref", help="generate reference triangle 'rate;x0,x1;y0'")
parser.add_argument("--group0", help="group input data by plot")
parser.add_argument("--comment", help="place some text in the bottom right corner, like the git hash, date, etc")
parser.add_argument("--marker", help="set the marker", default="+")
parser.add_argument("--figsize", help="figure size 'w,h' in inches", default="6,6")
parser.add_argument("--eoc", help="plot experimental order of convergence log(y_i/y_{i-1})/log(x_i/x_{i-1}) instead of y", action="store_true")
parser.add_argument("--tikz", help="prefix for pgfplots output: writes <prefix>.tex + <prefix>.csv; the .tex compiles standalone and is \\input-able via the standalone package; width via \\figurewidth")

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

for idx, name0, series, refs in collect_plots(df, args, tx, ty):
    if args.tikz: render_tikz(idx, name0, series, refs, args)
    render_matplotlib(idx, name0, series, refs, args, plot_func)

if not args.nshow:
    plt.show()
