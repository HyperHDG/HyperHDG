#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne12-01-conv-trace}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_trace, tau_s, e_abs}' $OUT.json \
    | experiments/plot.py -g deg,tau_s -x nx --trans '1/x,y' -y e_trace --eoc -w "deg<=2"\
                          --log x  --xbase 2 --xlabel 'discretization size h' \
                          --ylabel 'EOC of maximum $L^2$ error in trace variables'
#--ref "2;64,32;4e-3|3;64,32;3e-5|4;32,16;3e-6" \
#--save "$OUT.png,$OUT.pgf" --figsize '2.5,2.5'
#--legend ''

