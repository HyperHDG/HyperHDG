#!/bin/bash
# Follow-up to ne9-09-conv-xt.sh: the combined x+t run shows deg 5 (s = 3) limited to total
# order ~4 instead of 6. Separate the components: x-only arm (nt 256 fixed, temporal error
# negligible) and t-only arm (nx 16 fixed, spatial floor ~2e-8). First result: spatial is
# clean h^6+; the temporal arm carries the limit (~order 4, erratic), in both e_rel and
# e_trace -- rerun after the s = 3 endpoint-trace static solve lands.
# Rates: experiments/ne9-09-conv-xt-rates.sh output/ne9-09-conv-xt-01.json (nx column is
# meaningful for the x-arm; read the t-arm rows by nt).
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/complex/experiments
export OMP_NUM_THREADS=1

cmake --build --preset complex --target timowave

parallel --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg 5 -nx {1} -test wave4 -nt {2}" \
  ::: "2 256" "4 256" "8 256" "16 2" "16 4" "16 8" "16 16" "16 32"
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
