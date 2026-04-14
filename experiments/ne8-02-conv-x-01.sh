#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne8-02-conv-x}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, e_abs}' $OUT.json \
    | experiments/plot.py -x nx --trans '1/x,y' -y e_abs --log xy --save "$OUT.png" \
        --xlabel 'discretization size h' --ylabel 'maximum $L^2$ error'
