#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-04-procs-big}" # set default if unset
jq -c '.Stdout | .mpi.sz as $p | .net2as.load_bal | to_entries[] | {p: $p, load_bal: .value, type: .key}' $OUT.json \
  | experiments/plot.py -x p -y load_bal -g type --save $OUT.png,$OUT.pgf --log x
