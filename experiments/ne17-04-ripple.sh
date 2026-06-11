#!/bin/bash
set -xeo pipefail

# "ripple in the pond": a brief, localized tap at the center of a large clamped
# grid. Unlike ne17-03 (broad, slow tap -> the center stays excited and standing
# modes form), here the pulse is narrow in space and short in time relative to
# how far the wave travels, so the excitation leaves the center as an expanding
# ring and the center settles back to rest once the front has passed.
#
# levers (vs drumhead):
#   --grid 100   more cells between center and rim (domain is always unit square,
#                so "bigger pond" = finer grid + narrower pulse, not a bigger box)
#   std_x 3e-2   point-like splash (~3h, h=1/99); not a broad push
#   std_t 1.5e-2 brief impulse, over by t~0.1
#   tap_time 5e-2 (~3.3 std_t, so the tap ramps up from ~0 cleanly)
#   T 0.6        ~one center->rim crossing (c~1), stop before reflections return
# energy is fixed (amplitude = energy/(std_x^2 std_t)), so narrowing keeps the
# total impulse constant.

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
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py --grid 100 -o $DOMAIN --dirichlet xmin=63 xmax=63 ymin=63 ymax=63
$BUILD/timowave -domain $DOMAIN -test drumhead \
                -std_x 3e-2 -std_t 1e-2 -energy 5e-4 -tap_time 5e-2 -nt 75 -T 0.6 \
                -pc_type lu -ksp_type preonly \
                -plot $WAVE -deg 2 -log_view :$PERF:ascii_flamegraph \
                -print_timestep | tee $LOG
experiments/netvis.py $WAVE --view iso --beams 1 -o $AVI --color-by values:8
ffmpeg -i $AVI -c:v libx264 -c:a aac $VID
