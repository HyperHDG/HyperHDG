#!/bin/bash
# Temporal convergence of the Gauss stepping (s=1 implicit midpoint), replacing the retired
# trapezoidal theta scheme. Compare against the theta=0.5 rows of output/GOLDEN-ne9-01-conv-t.json:
# same order-2 rates and iteration counts; e_rel differs at O(dt^2) only (midpoint vs
# endpoint-averaged load sampling of the time-dependent Dirichlet data). Small arms on purpose.
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/openblas/experiments
export OMP_NUM_THREADS=1

cmake --build --preset openblas --target timowave

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave $DOMAIN -deg 6 -nx 20 -test wave4 -nt {1} -pc_type none" \
  ::: 8 16 32 64 128
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
