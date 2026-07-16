#!/bin/bash
# Temporal convergence sweep of the Gauss stepping. timowave.cxx couples the stage count to
# the spatial degree (deg 1/3/5 -> s = 1/2/3), so the matched design shows e_rel dropping at
# order dt^{2s} = 2/4/6 per nt doubling. nx is fixed per degree (256/64/32) so the spatial
# floor stays well below the temporal error at nt 32 (constants from the GOLDEN-ne9 tables).
# Status: deg 1 clean 2; deg 3 averages 4 with non-monotone steps (3.0..4.4 -- max-over-time
# e_rel artifact, identical in the old ne9-08 data); deg 5 (s = 3) sits at exactly order 4,
# pending the endpoint-trace static solve (spatial h^6 is clean per ne9-09-conv-x, and the
# limit shows in e_rel and e_trace alike). Stage systems are dense complex LU.
# Rates: experiments/ne9-conv-rates.sh output/ne9-10-conv-t.json nt
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
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -nt {3}" \
  ::: 1 3 5 :::+ 256 64 32 ::: 2 4 8 16 32
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
