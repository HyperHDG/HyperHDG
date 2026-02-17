#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=/tmp/np6630/HyperHDG

parallel -j 1 --progress --bar --results $OUT.json \
  "mpirun -n {1} build/rel/experiments/network -plot 0 -domain $DATA_DIR/domains/paper-big.geo.h5 -net2as_p 16 -net2as_pc_factor_mat_solver_type mumps -net2as_cb_type pu -mat_partitioning_type parhip -parhip_mode 1" \
    ::: 1 2 4 8 16 32 64 128
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
