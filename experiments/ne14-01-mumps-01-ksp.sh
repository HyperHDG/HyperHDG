#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-01-mumps}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {t_ksp, n: .net2as.global.size}' $OUT/res2.json \
    | experiments/plot.py -x n -y t_ksp --log xy --trans 6*x,y \
        --save "$OUT/ne14-01-mumps-01-ksp.png" \
        --ref "1.5;2e5,7e5;2e0" \
        --xlabel 'number of DOFs' --ylabel 'matrix factorization time'
#        --figsize '2.5,2.5' \
#        --legend ''

