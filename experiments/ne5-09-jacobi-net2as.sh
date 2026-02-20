#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="paper-small"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN.geo.h5 -mat_cache $DATA_DIR/output/$DOMAIN.bin"
BIN_DIR=build/rel/experiments
MUMPS="-net2as_pc_factor_mat_solver_type mumps"

export OMP_NUM_THREADS=1

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/network $DOMAIN -ksp_monitor_yaml -pc_type {1} -net2as_p {2}" \
  ::: jacobi net2as net2as :::+ "" 1 8
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}

