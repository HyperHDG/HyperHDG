#!/bin/bash
# Spatial convergence sweep of the Gauss stepping. timowave.cxx couples the stage count to
# the spatial degree (deg 1/3/5 -> s = 1/2/3), so the matched design shows e_rel dropping at
# order h^{p+1} = 2/4/6 per nx doubling. nt 256 is fixed: the temporal error dt^{2s} sits
# 100x+ below the spatial one at every level for all three degrees (constants from the
# GOLDEN-ne9 tables and ne9-09-conv-xt-01). Stage systems are dense complex LU.
# Rates: experiments/ne9-conv-rates.sh output/ne9-09-conv-x.json nx
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/complex/experiments
export OMP_NUM_THREADS=1

cmake --build --preset complex --target timowave

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -nt 256" \
  ::: 1 3 5 ::: 2 4 8 16 32
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
