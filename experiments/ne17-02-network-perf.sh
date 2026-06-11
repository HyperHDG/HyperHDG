#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
TRACE=$OUT/trace.h5
STATIC=$OUT/static.vtkhdf
WAVE=$OUT/wave.vtkhdf
AVI=$OUT/$NAME.avi
VID=$OUT/$NAME.mp4
IMG=$OUT/$NAME.png
PERF=$OUT/perf.flamegraph
LOG=$OUT/log.yaml
CUT=.5
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network --target timowave
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py -i $INPUT --clamp-xy $CUT -o $DOMAIN --dirichlet xmin=63 xmax=68
mpirun -n 2 $BUILD/network -domain $DOMAIN -strain .15 -pc_type cholesky -ksp_type preonly -pc_factor_mat_solver_type mumps \
                -plot $STATIC -log_view :$PERF:ascii_flamegraph \
                | tee $LOG.1
$BUILD/timowave -domain $DOMAIN -domain $DOMAIN -mat_only -test stiffness | tee $LOG.2
experiments/netvis.py $STATIC --view iso --beams 1 -o $IMG

