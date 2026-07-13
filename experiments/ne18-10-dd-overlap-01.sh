#!/bin/bash
set -xeo pipefail

: ${OUTDIR:=output}
OUT=$OUTDIR/ne18-10-dd-overlap
LOG=$OUT/data.json
IMG=$OUT/ne18-01-dd-10-overlap-01.png

jq -c '.Stdout | {
  label, iters: .iterations, cond, t_ksp, t_iteration,
  n_max: [.net2as.local[].size]|max,
  n_min: [.net2as.local[].size]|min,
  n_sum: [.net2as.local[].size]|add,
}' $LOG
#  | experiments/plot.py -x p -y cond -g nocoarse --save $IMG --log xy
