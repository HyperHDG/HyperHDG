#!/bin/bash
# Degree-parity probe for the trace superconvergence excess: deg 3/5 measure e_trace EOCs of
# ~5.7 / 7.5 (> p+2) on BOTH the compatible sin-profile and the stiff cos-profile arms
# (ne9-12-conv-x-hom-03, ne9-13-conv-trace-stiff), while deg 1 sits at p+2 = 3 and the old
# ne12-02 evidence for clean p+2 rests on p = 2. Missing data point: an EVEN degree in the
# new code. deg 2 with -stages 2 (order-4 time, compiled combo (2,2)) on the compatible arm,
# nt 1024 (temporal floor ~1e-11, trace at nx 128 ~1e-8 stays clear).
# Expect: clean 4.0 = p+2 -> parity pattern (even p exact, odd p gains ~1/2..1 order);
# ~4.5+ -> the excess is generic in the new apparatus and the old p=3 data was floor-warped.
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

parallel -j 6 --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave $DOMAIN -deg 2 -stages 2 -nx {1} -test wave4 -wave4_px $PX -nt 1024" \
  ::: 8 16 32 64 128
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
