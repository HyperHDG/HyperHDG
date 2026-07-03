#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/rel/experiments
export OMP_NUM_THREADS=1

cmake --build --preset rel --target timowave

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave $DOMAIN -deg {3} -nx {1} -theta .5 -test 4 -nt {4} -pc_type none -tau_s {2}" \
  ::: 2 4 8 16 32 64 ::: 1 0 -1 ::: 1 2 3 :::+ 500 1000 8000
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
