#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-03-conv-tau}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, e_abs, tau_s}' $OUT.json \
    | experiments/plot.py -x nx --trans '1/x,y' -y e_abs --log xy -g tau_s \
        --xbase 2 --xlabel 'discretization size h' --ylabel 'maximum $L^2$ error' \
        --figsize '2.5,2.5' --save "$OUT.png,$OUT.pgf" \
        --ref "2;64,32;6e-3|1;64,32;7e-2|2;32,16;1e-3" --legen ''
#--ref "2;64,32;4e-3|3;64,32;3e-5|4;32,16;3e-6"
# --legend ''
