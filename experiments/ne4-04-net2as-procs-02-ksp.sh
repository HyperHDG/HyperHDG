#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-04-net2as-procs}" # set default if unset
mkdir -p $OUT_DIR

jq -c '.Stdout | {"p": .mpi.sz, t_ksp}' $OUT.json \
    | jq -sc '.[0].t_ksp as $t1 | .[] | {p, E_ksp: ($t1 / (.p * .t_ksp))}' \
    | experiments/plot.py -f json -x p -y E_ksp --log x --save $OUT.png,$OUT.pgf
