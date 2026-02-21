#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-03-p-big}" # set default if unset
jq -c '.Stdout
  | .net2as.sz as $p | .net2as.cb_type as $cb | .iterations as $it
  | {t_ksp, t_it: .t_iteration, t: (.t_ksp + .t_iteration)}
  | to_entries[]
  | {p: $p, t: .value, type: .key, it: $it, cb: $cb}
' $OUT.json \
  | experiments/plot.py -x p -y t -g type,cb --log xy
