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
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
PERF=$OUT/perf.flamegraph
LOG=$OUT/log.json
DIRECT="-pc_type cholesky -ksp_type preonly -pc_factor_mat_solver_type mumps"
NET="-pc_type net2as -net2as_cb_type q1 -net2as_pc_factor_mat_solver_type mumps"
CUT=.4
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py --grid 100 -o $DOMAIN --dirichlet xmax=68 xmin=63
for nc in "" "-net2as_nocoarse"; do
for p in 1 2 4 8 16; do
  mpirun -n 8 $BUILD/network -domain $DOMAIN  -comp 2 -strain .15 $NET -net2as_p $p $nc \
                 -trace_view hdf5:$TRACE | yq -o json -I0 >> $LOG
done
done
