#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne8-01-conv-t}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {dt, e_abs}' $OUT.json \
    | experiments/plot.py -x dt -y e_abs --log xy --ref "2;1e-2,2e-2;1e-3" --save "$OUT.png,$OUT.pgf" \
        --xlabel 'timestep size $\Delta t$' --ylabel 'maximum $L^2$ error'
