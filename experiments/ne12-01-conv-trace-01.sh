#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne12-01-conv-trace}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_trace, tau_s, e_abs}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_trace -w 'tau_s==0'\
        --log xy  --xbase 2 --xlabel 'discretization size h' \
        --ylabel 'max $L^2$ error in $\lambda$' \
        --ref "2.5;64,32;2e-3|3.5;16,8;2e-4|4.5;16,8;1.5e-6" \
        --save "$OUT.png,$OUT.pgf" --tikz "$OUT" --legend ''
