#!/bin/bash
set -xeo

BUILD=build/rel/experiments
N=7
OUT=output/0529-norescale
INPUT=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/
DOMAIN=$OUT/domain.geo.h5
TRACE=$OUT/trace.h5
STATIC=$OUT/static.vtkhdf
WAVE=$OUT/wave.vtkhdf
VID=$OUT/0529-norescale.avi
IMG=$OUT/0529-norescale.png
NET="-pc_type net2as -net2as_p 1 -net2as_cb_type pu -net2as_pc_factor_mat_solver_type mumps"
export OMP_NUM_THREADS=1

mkdir -p $OUT
cmake --build --preset rel
python experiments/make_geo2.py -i $INPUT --clamp-xy .2 -o $DOMAIN --dirichlet xmax=68 xmin=63
$BUILD/network -domain $DOMAIN  -comp 2 -strain .15 $NET \
               -trace_view hdf5:$TRACE -plot $STATIC
$BUILD/timowave -domain $DOMAIN -comp 2 -strain .15 -nt 100 -T 1e-5 $NET \
                -static $TRACE -plot $WAVE -deg 3 -net2as_wave -tau 1e0
experiments/netvis.py $WAVE --view iso -o $VID -r 1
experiments/netvis.py $WAVE --view pside --frames 0,.5e-5,1e-5 --frame-colors tan,purple,blue \
                      -o $IMG --axis 0 --beams 1
#ffmpeg -i $OUT/$NAME-wave.avi -c:v libx264 -c:a aac $OUT/$NAME-wave.mp4
cp $OUT $OUT.$(date +%s)
