#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-05-net2as-subdoms}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout | .net2as.sz as $p | .net2as.local[] | {p: $p, fill, time, size}' $OUT.json \
    | experiments/plot.py -f json -x size -y fill -g p --log xy --ref ".3;2e4,6e4;2"

