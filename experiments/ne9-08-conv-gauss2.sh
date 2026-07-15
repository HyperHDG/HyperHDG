#!/bin/bash
# Convergence of the 2-stage Gauss stepping (deg 3 -> s = 2 automatically; order 4 in time).
# Temporal arm: deg 3, nx 32, nt 4..64 -- e_rel ratios ~16 per nt-doubling until the spatial
# floor (~3.07e-6, the GOLDEN-ne9-02 deg-3 nx-32 value). Spatial arm: nt 64 suffices for h^4
# down to nx 16 (the golden CN table needed nt 8000); nx 32 runs at nt 128 to clear the floor.
# Stage systems are dense complex LU (no PETSc solver options needed).
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/openblas/experiments
export OMP_NUM_THREADS=1

cmake --build --preset openblas --target timowave

parallel --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg 3 -nx {1} -test wave4 -nt {2}" \
  ::: "32 4" "32 8" "32 16" "32 32" "32 64" "2 64" "4 64" "8 64" "16 64" "32 128"
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
