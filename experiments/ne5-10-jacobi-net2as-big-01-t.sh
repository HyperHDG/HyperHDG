#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-10-jacobi-net2as-big}" # set default if unset
jq -c '.Stdout
  | .net2as.sz as $p
  | {"net2as": "cg+as", "jacobi": "cg+jacobi"}[.pc_type] as $type
  | if .net2as.sz == 1 then "direct" else $type end as $type
  | .ksp_monitor[]
  | {p: $p // 1, type: $type, time, rnorm}
' $OUT.json \
  | experiments/plot.py -x time -y rnorm -g type --log y --marker '' --save $OUT.png,$OUT.pgf --ylabel 'rel. residual'
