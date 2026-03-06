#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne6-TestTimoWave4-p2-nt}" # set default if unset
mkdir -p $OUT_DIR

jq -cs 'map(.Stdout | {e_abs, nx, nt}) | sort_by(.nt).[]' $OUT.json
