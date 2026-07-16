#!/bin/bash
# Trace-superconvergence probe for the compatible (sin-profile) wave4 arm: the dyadic sweeps
# measure e_trace EOCs ABOVE the expected p+2 (deg 3: ~5.8, deg 5: ~7.5), unlike the old
# cos-profile data (clean p+2). Two diagnostics:
#  (1) ternary refinement nx = 9/27/81 (deg 3): interior nodes never coincide with the zeros
#      of sin(2 pi s) -- if the excess rate is a node/zero alignment artifact of dyadic
#      refinement, p+2 = 5 reappears here (EOC = log3 ratios!);
#  (2) deep deg-1 points nx = 64..256 at nt 4096 (s = 1, temporal floor ~8e-7): the shallow
#      sweep's deg-1 trace EOCs were too noisy to read (2.2 .. 3.7).
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

parallel -j 6 --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -wave4_px $PX -nt {3}" \
  ::: "3 9 1024" "3 27 1024" "3 81 1024" "1 64 4096" "1 128 4096" "1 256 4096"
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
