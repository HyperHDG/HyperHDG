#!/bin/bash
# Simultaneous x+t convergence of the Gauss stepping: nx and nt are refined together (dt ~ h).
# timowave.cxx couples the stage count to the spatial degree (deg 1/3/5 -> s = 1/2/3), so the
# spatial h^{p+1} and temporal dt^{2s} components refine at the same rate and e_rel drops with
# total order 2 / 4 / 6 per doubling (deg 5 is pre-asymptotic at nx 2; clean from nx 4 on).
# Stage systems are dense complex LU (no PETSc solver options needed).
# Rates table: experiments/ne9-09-conv-xt-rates.sh output/ne9-09-conv-xt.json
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
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -nt {2}" \
  ::: 1 3 5 ::: 2 4 8 16 32
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
