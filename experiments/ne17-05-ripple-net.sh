#!/bin/bash
set -xeo pipefail

# "ripple in the pond" on the disordered fibre network, in the SAME regime as
# ne17-06-sinclamp-net.sh (MHz drive, micron domain, microsecond run). A brief,
# localized out-of-plane tap at the centre of the morgan network, clamped on all
# four borders. The pulse is narrow in space and short in time relative to how
# far the wave travels, so it leaves the centre as an expanding (and strongly
# scattered) ring rather than building up standing modes.
#
# regime (mirrors ne17-06):
#   --clamp .4         same clamped sub-network (x-extent ~3.2 mm)
#   std_x 1e2 (100um)  point-like splash: ~1 flexural wavelength, ~3% of domain
#   std_t 5e-8 (50ns)  pulse bandwidth ~1/(2*pi*std_t) ~ 3 MHz  (== ne17-06 -freq 3e6)
#   tap_time 1.5e-7    ~3*std_t, so the tap ramps up from ~0 cleanly
#   T 1e-6, nt 60      ~one centre->rim crossing (c~1.7e9 um/s), == ne17-06
# energy is fixed (amplitude = energy/(std_x^2 * std_t)); the solver is linear, so
# -energy only rescales the result -- tune it for a visible warp amplitude.

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
CUT=.2
MUMPS="-pc_type cholesky -ksp_type preonly -pc_factor_mat_solver_type mumps"
NET="-pc_type net2as -net2as_p 1 -net2as_cb_type pu -net2as_pc_factor_mat_solver_type cholmod"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py -i $INPUT --clamp $CUT -o $DOMAIN --dirichlet xmin=63 xmax=63 ymin=63 ymax=63 -t 1e-1
mpirun -n 8 $BUILD/timowave -domain $DOMAIN -test drumhead \
                -std_x 3e1 -std_t 8e-8 -energy 3e1 -tap_time 2e-7 -T 2e-6 -nt 60 $NET \
                -plot $WAVE -deg 6 -log_view :$PERF:ascii_flamegraph -tau 1e0 \
                -print_timestep | tee $LOG
experiments/nethist.py $WAVE
#experiments/netvis.py $WAVE --view iso --beams 1 -o $AVI --color-by values:8
#ffmpeg -i $AVI -c:v libx264 -c:a aac $VID
#experiments/nethist.py $WAVE --bins 100
