#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-09-jacobi-net2as}" # set default if unset
jq -c '.Stdout | {p: .net2as.sz, t: (.t_ksp + .t_iteration), type: if .net2as.sz == 1 then "direct" else .pc_type end}' $OUT.json
