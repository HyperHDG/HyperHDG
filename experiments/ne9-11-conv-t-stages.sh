#!/bin/bash
# Temporal order of the Gauss stepping with the SPATIAL DISCRETIZATION HELD CONSTANT
# (deg 5, nx 256, cross2): only the stage count varies (-stages 1/2/3, decoupled from the
# degree). Two arms:
#   stiff:       wave4 as-is -- time-dependent Dirichlet data drives every boundary node, the
#                delta-like spatial coupling excites the stiff modes (dt*|lambda| >> 1), and
#                Gauss collocation reduces to its stiff order: ~s+1 (odd s) / s (even s).
#                Expect ~2 / 2..3 / 4 -- NOT 2s, independent of the spatial resolution.
#   homogeneous: wave4 with -wave4_px -pi/2 -- the spatial standing-wave factor becomes
#                sin(w s), zero at all tips: Dirichlet data identically zero, only the
#                spatially smooth interior forcing remains, no stiff excitation. Expect the
#                classical orders 2s = 2 / 4 / 6 (deg5-nx256 spatial floor ~1e-13).
# Mechanism isolated in experiments/gauss_value_form_check.py (Prothero-Robinson arm).
# Rates: experiments/ne9-conv-rates.sh output/ne9-11-conv-t-stages-{stiff,hom}.json nt gauss_stages
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/complex/experiments
export OMP_NUM_THREADS=1
PX=-1.5707963267948966

cmake --build --preset complex --target timowave

# -j 4: the deg-5 nx-256 per-edge complex LU caches are ~0.5 GB per job
parallel -j 4 --progress --bar --results $OUT-stiff.json \
  "$BIN_DIR/timowave $DOMAIN -deg 5 -nx 256 -stages {1} -test wave4 -nt {2}" \
  ::: 1 2 3 ::: 4 8 16 32 64
echo "gen exit (stiff): $?"
parallel -j 4 --progress --bar --results $OUT-hom.json \
  "$BIN_DIR/timowave $DOMAIN -deg 5 -nx 256 -stages {1} -test wave4 -wave4_px $PX -nt {2}" \
  ::: 1 2 3 ::: 4 8 16 32 64
echo "gen exit (hom): $?"
yq -i '.Stdout |= from_yaml' $OUT-stiff.json
yq -i '.Stdout |= from_yaml' $OUT-hom.json
cp $OUT-stiff.json{,.$(date +%s)}
cp $OUT-hom.json{,.$(date +%s)}
