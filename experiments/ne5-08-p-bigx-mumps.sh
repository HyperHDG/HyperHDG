#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=/tmp/np6630/HyperHDG
#DOMAIN="paper-big"
#DOMAIN="-domain $DATA_DIR/domains/$DOMAIN.geo.h5 -mat_cache $DATA_DIR/output/$DOMAIN.bin"
BIN_DIR=build/rel/experiments

mkdir -p $DATA_DIR
mkdir $DATA_DIR/{output,domains}

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/network -domain $DATA_DIR/domains/paper-big{1}.geo.h5 -mat_cache $DATA_DIR/output/paper-big{1}.mat.bin -net2as_p {2} -net2as_cb_type {3} -net2as_pc_factor_mat_solver_type mumps -net2as_print_local" \
  ::: 2 3 4 5 ::: 1 2 4 8 16 24 ::: pu
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}

