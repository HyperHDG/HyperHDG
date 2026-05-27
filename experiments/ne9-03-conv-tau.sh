#!/bin/bash
set -xeu pipefail

: ${OUT_DIR:=output}
: ${OUT:=$OUT_DIR/$(basename ${0%.sh})}
mkdir -p $OUT_DIR

DOMAIN="cross2.geo"
export DOMAIN="-domain domains/$DOMAIN"
export BIN_DIR=build/rel/experiments
export OMP_NUM_THREADS=1

parallel --progress --bar --results $OUT.json \
  '$BIN_DIR/timowave $DOMAIN -deg {1} -nt {2} -nx {3} -tau {=3 $_ = 2**-$_ =} -theta .5 -test 4 -pc_type none' \
  ::: 1 :::+ 500 ::: 2 4 8 16 32 64
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
