#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-03-threaded}"
mkdir -p $OUT_DIR

jq -cs '
  map(.Stdout | {t_ksp, threads, n: .net2as.global.size, type: .net2as.local[0].factor_type})
  | (INDEX(.[] | select(.threads==1); .type)) as $b
  | .[] | .t_ksp = ($b[.type].t_ksp / .t_ksp) / .threads
' "$OUT/res2.json" \
    | experiments/plot.py -x threads -y t_ksp -g type --log x \
        --save "$OUT/ne14-03-threaded-02-ef.png" \
        --xlabel 'number of threads' --ylabel 'parallel efficiency' \
        --xbase 2
