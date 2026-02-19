#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-04-procs-big}" # set default if unset
jq -sc '
  map(.Stdout | {p: .mpi.sz, t_t2f, t_ksp, t_it: .t_iteration}) | sort_by(.p) |
  . as $all |
  ["t_t2f", "t_ksp", "t_it"] | .[] as $key |
  $all | (.[0][$key]) as $base |
  .[] | {p, stage: $key, E: ($base / (.[$key] * .p))}
' $OUT.json | experiments/plot.py -x p -y E -g stage --log x
