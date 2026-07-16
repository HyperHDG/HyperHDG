#!/bin/bash
# Stabilization plot of ne9-03-conv-tau.sh: e_rel vs h at p = 3, one line per tau ~ h^s
# exponent s in {-1, 0, 1}. Right panel of the manuscript's shared-axis row: y limits pinned,
# tick labels suppressed (they live on the left neighbor). Writes $OUT.png (view: sxiv) and
# the standalone pgfplots pair $OUT.tex/$OUT.csv.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-03-conv-tau}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, s: .tau_s, e_rel}' $OUT.json \
    | experiments/plot.py -g s -x nx --trans '1/x,y' -y e_rel \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylim '1e-11,10' --ytickoff \
        --ref "3;32,16;2e-5|4;32,16;8e-7|5;32,16;2e-8" \
        --save "$OUT.png" --tikz "$OUT" --legend 'lower right' --nshow
