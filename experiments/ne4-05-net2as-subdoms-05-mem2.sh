#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-05-net2as-subdoms}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout.net2as | {p: .sz, nz: (.local | map(.nz_fac) | add)}' $OUT.json \
  | experiments/plot.py -f json -x p -y nz --log xy
