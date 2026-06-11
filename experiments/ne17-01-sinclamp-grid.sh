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
CUT=1
MUMPS="-pc_type cholesky -ksp_type preonly -pc_factor_mat_solver_type mumps"
NET="-pc_type net2as -net2as_p 1 -net2as_cb_type pu -net2as_pc_factor_mat_solver_type cholmod"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py --grid 20 -o $DOMAIN --dirichlet xmin=63
mpirun -n 2 $BUILD/timowave -domain $DOMAIN -strain .03 -freq 3 -nt 30 $NET \
                -plot $WAVE -deg 2 -log_view :$PERF:ascii_flamegraph \
                -print_timestep -test sinclamp | tee $LOG
experiments/netvis.py $WAVE --view iso --beams 1 -o $AVI
ffmpeg -i $AVI -c:v libx264 -c:a aac $VID
