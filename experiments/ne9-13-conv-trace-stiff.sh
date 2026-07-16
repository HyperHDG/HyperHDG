#!/bin/bash
# Trace convergence on the ORIGINAL (stiff, cos-profile) wave4: does the hybrid variable
# recover the clean p+2 superconvergence of the old ne12-02 data, in contrast to the
# compatible sin-profile arm (~p+2.5..3, see ne9-12-conv-x-hom-03)? If yes, the excess is a
# parity/profile cancellation, not a property of the method.
# Time stepping on the stiff arm is order-reduced (s+1 / s), so nt is chosen per degree to
# push the temporal e_trace floor below the finest spatial trace error:
#   deg 1 (s=1, order 2):     nt 8192  -> floor ~2e-9  vs trace(nx 256) ~ 6e-6
#   deg 3 (s=2, order ~2.5):  nt 16384 -> floor ~1e-11 vs trace(nx 64)  ~ 3e-10
#   deg 5 (s=3, order 4):     nt 2048  -> floor ~5e-14 vs trace(nx 16)  ~ 1e-10 (roundoff
#                             ~5e-12 caps nx anyway)
# e_rel is temporally floored on this arm and NOT usable from this run; only e_trace is.
# Rates: yq per-degree ratios, or ne9-conv-rates.sh (nx doublings within each degree).
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/complex/experiments
export OMP_NUM_THREADS=1

cmake --build --preset complex --target timowave

joblist() {
  for nx in 16 32 64 128 256; do echo "1 $nx 8192";  done
  for nx in 4 8 16 32 64;     do echo "3 $nx 16384"; done
  for nx in 2 4 8 16;         do echo "5 $nx 2048";  done
}

parallel -j 6 --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -nt {3}" \
  :::: <(joblist)
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
