#!/bin/bash
# Stabilization sweep (recreated for the new Gauss stepping): tau ~ h^s for s in {1, 0, -1}
# (-tau_s) at fixed p = 3 (matched s = 2, temporal order 4) on the COMPATIBLE wave4 arm
# (-wave4_px -pi/2, homogeneous Dirichlet -- no stiff order reduction, see
# gauss_order_reduction.tex). nt 1024 keeps the temporal floor ~1e-11, below the spatial
# error through nx 256 (only the tau~1/h nx-256 point touches it). Successor of the CN-era
# ne9-03-conv-tau (theta scheme, stiff wave4), whose per-(deg,nx) nt gymnastics the order-4
# stepping makes unnecessary.
# Measured orders: tau~h -> 3 (= p), tau=1 -> 4 (= p+1); tau~1/h ALSO tends to 4: its EOCs
# decline 5.1 -> 4.4 over nx 16..128 -- an h^5 pre-asymptotic transient on top of the
# regular h^4 term, NOT superconvergence.
# Rates: experiments/ne9-conv-rates.sh output/ne9-03-conv-tau.json nx tau_s
# Plot:  experiments/ne9-03-conv-tau-01.sh
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
  "$BIN_DIR/timowave $DOMAIN -deg 3 -nx {1} -tau_s {2} -test wave4 -wave4_px $PX -nt 1024" \
  ::: 2 4 8 16 32 64 128 256 ::: 1 0 -1
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
