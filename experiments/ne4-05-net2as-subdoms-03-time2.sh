#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-05-net2as-subdoms}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout | {p: .net2as.sz, t_ksp, t_iteration}' $OUT.json

