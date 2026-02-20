#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="paper-big"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN.geo.h5"
BIN_DIR=build/rel/experiments
NET="-net2as_p 24 -net2as_print_local"
FACTYPE="-net2as_pc_factor_mat_solver_type" # only affects the local mats

parallel -j 1 --progress --bar --results $OUT.json \
  "mpirun -n {1} $BIN_DIR/network -plot 0 $DOMAIN $NET -net2as_cb_type {2} $FACTYPE {3}" \
  ::: 1 2 4 8 16 32 64 128 ::: q1 pu ::: petsc mumps
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}

