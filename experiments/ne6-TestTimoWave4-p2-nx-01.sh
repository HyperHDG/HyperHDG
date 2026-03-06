#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne6-TestTimoWave4-p2-nx}" # set default if unset
mkdir -p $OUT_DIR

jq -cs 'map(.Stdout | {e_abs, nx, nt}) | sort_by(.nx).[]' $OUT.json \
   | experiments/plot.py -x nx -y e_abs --log xy --xbase 2 --ref '-3;8,32;6e-2'
