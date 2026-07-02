#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-01-conv-t}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {dt, e_rel, e_abs, theta}' $OUT.json \
    | experiments/plot.py -g theta -x dt -y e_rel --log xy \
        --xbase 2 --figsize '2.5,2.5' \
        --ref "2;.015625,.03125;2e-3|1;1e-2,2e-2;8e-2" \
        --save "$OUT.png,$OUT.pgf" --tikz "$OUT" \
        --xlabel 'timestep size $\Delta t$' --ylabel 'maximum $L^2$ error' \
        --legend ''
