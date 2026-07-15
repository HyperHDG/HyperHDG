#!/bin/bash
# Spatial convergence h^{p+1} of the Gauss stepping (s=1 implicit midpoint). Compare against
# output/GOLDEN-ne9-02-conv-x.json (theta=0.5): same rates; small deviations from the O(dt^2)
# load sampling difference. deg 3 / nt 8000 arm of the golden omitted to keep the run short.
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
  "$BIN_DIR/timowave $DOMAIN -deg {2} -nx {1} -test wave4 -nt {3} -ksp_type preonly -pc_type lu" \
  ::: 2 4 8 16 32 ::: 1 2 :::+ 500 1000
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
