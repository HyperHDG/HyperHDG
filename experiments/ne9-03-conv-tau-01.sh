#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-03-conv-tau}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, e_rel, tau_s}' $OUT.json \
    | experiments/plot.py -x nx --trans '1/x,y' -y e_rel --log xy -g tau_s \
        --xbase 2 --xlabel 'discretization size h' --ylabel 'maximum $L^2$ error' \
        --figsize '2.5,2.5' --save "$OUT.png,$OUT.pgf" --tikz "$OUT" \
        --ref "2;64,32;4e-3|1;64,32;9e-2|2;32,16;1e-3" --legend ''
