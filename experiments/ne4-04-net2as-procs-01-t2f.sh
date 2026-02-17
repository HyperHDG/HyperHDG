#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-04-net2as-procs}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout | {"p": .mpi.sz, t_t2f}' $OUT.json \
    | jq -sc '.[0].t_t2f as $t1 | .[] | {p, E_t2f: ($t1 / (.p * .t_t2f))}' \
    | experiments/plot.py -f json -x p -y E_t2f --log x
