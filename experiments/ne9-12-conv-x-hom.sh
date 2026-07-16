#!/bin/bash
# Spatial convergence sweep on the COMPATIBLE (homogeneous-Dirichlet) wave4 arm:
# -wave4_px -pi/2 puts the standing-wave nodes on every tip, so the Dirichlet data is
# identically zero and the stepping delivers its full temporal order 2s (no stiff order
# reduction, see gauss_order_reduction.tex). Matched stages (deg 1/3/5 -> s = 1/2/3),
# nt 256 fixed: the temporal error sits far below the spatial one at every level.
# Expect clean h^{p+1} = 2/4/6.
# Rates: experiments/ne9-conv-rates.sh output/ne9-12-conv-x-hom.json nx
# Plot:  experiments/ne9-12-conv-x-hom-01.sh
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

parallel --progress --bar --results $OUT.json \
  "$BIN_DIR/timowave $DOMAIN -deg {1} -nx {2} -test wave4 -wave4_px $PX -nt 256" \
  ::: 1 3 5 ::: 2 4 8 16 32
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
