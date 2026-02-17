#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne4-04-net2as-procs}" # set default if unset
mkdir -p $OUT_DIR

yq -I0 -o=json '.Stdout | {"p": .mpi.sz, "t_t2f": .t_t2f}' $OUT.json | experiments/plot.py -f json -x p -y t_t2f
