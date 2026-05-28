#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne12-01-conv-trace}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_trace, tau_s, e_abs}' $OUT.json \
    | experiments/plot.py -g deg,tau_s -x nx --trans '1/x,y' -y e_trace \
        --log xy  --xbase 2 --xlabel 'discretization size h' \
        --ylabel 'maximum $L^2$ error in trace variables' \
        --ref "2.5;64,32;4e-3|3.5;64,32;3e-6|4.5;16,8;3e-6" \
        --save "$OUT.png"
