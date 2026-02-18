#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="paper-small"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN.geo.h5"
BIN_DIR=build/rel/experiments

parallel -j 1 --progress --bar --results $OUT.json \
  "mpirun -n {1} $BIN_DIR/network -plot 0 $DOMAIN -net2as_p 16 -net2as_print_local" \
  ::: 1 2 4 8 16 32 64
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}

