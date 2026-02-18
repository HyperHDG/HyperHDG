#!/bin/bash
set -x
: "${OUT:=${OUT_DIR:=output}/ne5-06-p-big-mumps}" # set default if unset
jq -c '.Stdout | {p: .net2as.sz, t_ksp, t_it: .t_iteration, t: (.t_ksp + .t_iteration), it: .iterations, cb: .net2as.cb_type}' $OUT.json
