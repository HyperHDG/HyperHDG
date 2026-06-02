#!/bin/bash
set -xeo pipefail

: ${BUILD:=build/rel/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
TRACE=$OUT/trace.h5
STATIC=$OUT/static.vtkhdf
WAVE=$OUT/wave.vtkhdf
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
PERF=$OUT/perf.flamegraph
LOG=$OUT/log.yaml
NET="-pc_type net2as -net2as_p 1 -net2as_cb_type pu -net2as_pc_factor_mat_solver_type cholmod"
CUT=.5
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset rel --target network
cp $BUILD/network experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo '-dirty' >> $OUT/rev; fi

python experiments/make_geo2.py -i $INPUT --clamp-xy $CUT -o $DOMAIN --dirichlet xmax=68 xmin=63
mpirun -n 2 $BUILD/network -domain $DOMAIN  -comp 2 -strain .15 $NET \
               -trace_view hdf5:$TRACE -plot $STATIC
experiments/netvis.py $STATIC --view iso --beams 1
#experiments/netvis.py $WAVE --view pside --frames 0,.5e-5,1e-5 --frame-colors tan,purple,blue \
#                      -o $IMG --axis 0 --beams 1
#ffmpeg -i $OUT/$NAME-wave.avi -c:v libx264 -c:a aac $OUT/$NAME-wave.mp4
