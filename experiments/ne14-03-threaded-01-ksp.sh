#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-03-threaded}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {t_ksp, threads, n: .net2as.global.size, type: .net2as.local[0].factor_type}' "$OUT/res2.json" \
    | experiments/plot.py -x threads -y t_ksp -g type \
        --save "$OUT/ne14-03-threaded-01-ksp.png" \
        --xlabel 'number of threads' --ylabel 'matrix factorization time'
#        --ref "1.5;2e5,7e5;2e0" \
#        --figsize '2.5,2.5' \
#        --legend ''

