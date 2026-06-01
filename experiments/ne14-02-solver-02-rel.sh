#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-02-solver}"
mkdir -p $OUT_DIR

jq -cs '
  map(.Stdout | {t_ksp, n: .net2as.global.size,
                 type: .net2as.local[0].factor_type, blas})
  | (INDEX(.[] | select(.blas=="openblas" and .type=="cholmod"); .n)) as $b
  | .[] | .t_ksp /= $b[.n|tostring].t_ksp
' "$OUT/res2.json" \
    | experiments/plot.py -x n -y t_ksp -g blas,type --log xy --trans 6*x,y \
        --save "$OUT/ne14-02-solver-01-ksp.png" \
        --xlabel 'number of DOFs' --ylabel 'matrix factorization time'
#        --ref "1.5;2e5,7e5;2e0" \
#        --figsize '2.5,2.5' \
#        --legend ''

