#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="line1.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/rel/experiments

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave -plot 0 $DOMAIN -deg 2 -nt {1} -theta .5 -test 4 -tau 1 -nx 500" \
  ::: 1 2 4 8 16 32 64 128 256 512
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
