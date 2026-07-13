#!/bin/bash
set -xeo pipefail

: ${OUTDIR:=output}
OUT=$OUTDIR/ne18-04-dd-cutoff
LOG=$OUT/log.json
IMG=$OUT/ne18-04-dd-cutoff-01-p-cond.png

jq -cs 'sort_by(.net2as.nocoarse, .net2as.p[0]) | .[]
       | {p: .net2as.p[0], nocoarse: .net2as.nocoarse, iters: .iterations, cond, creason}' $LOG \
  | experiments/plot.py -x p -y cond -g nocoarse --save $IMG --log xy --xbase 2
