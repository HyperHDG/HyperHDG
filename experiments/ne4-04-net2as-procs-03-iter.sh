#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-04-net2as-procs}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout | {"p": .mpi.sz, t_iteration}' $OUT.json \
    | jq -sc '.[0].t_iteration as $t1 | .[] | {p, E_iter: ($t1 / (.p * .t_iteration))}' \
    | experiments/plot.py -f json -x p -y E_iter --log x
