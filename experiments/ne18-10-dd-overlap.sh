#!/bin/bash
# Sweep the net2as "pu" overlap distance (as a fraction of the weighted graph diameter) at a fixed
# number of subdomains, plus the q1-trim baseline (same subdomain count). For every config we plot
# the residual history (-ksp_monitor_yaml) against wall-clock time. The monitor's time field already
# includes the factorization (it is logged at iteration 0), so configurations that trade a cheaper
# factorization for more iterations are compared fairly. q1-trim is included because it produces p*p
# subdomains, matching pu.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
: ${P:=4}                               # subdomains per axis -> p*p subdomains
: ${FRACS:="0.1 0.2 0.3"}      # overlap distance / per-subdomain weighted diameter
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
MAT=$OUT/mat.bin
LOG=$OUT/data.json
PLOT=$OUT/$NAME.png
CUT=.4
# unpreconditioned norm so every config converges in the same true residual ||b-Ax|| (the
# preconditioned norm differs per preconditioner, making the configs stop at different residuals)
NET="-pc_type net2as -net2as_pc_factor_mat_solver_type mumps -ksp_norm_type unpreconditioned"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/{make_geo2,plot}.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py -i $INPUT --clamp-xy $CUT -o $DOMAIN \
       --dirichlet xmax=63 xmin=63 ymin=63 ymax=63

cmd() {
  label=$1; shift
  echo "echo \"label: $label\"; $BUILD/network -mat_cache $MAT -test constant -domain $DOMAIN $NET \
         -net2as_p $P -net2as_print_local -ksp_monitor_yaml -mem_max $@"
}

parallel --progress --bar --results $LOG \
  <<EOF
$(cmd "q1"            -net2as_cb_type q1 -net2as_cb_trim)
$(cmd "pu delta=5e-2" -net2as_cb_type pu -net2as_overlap_frac 5e-2)
$(cmd "pu delta=1e-1" -net2as_cb_type pu -net2as_overlap_frac 1e-1)
$(cmd "pu delta=2e-1" -net2as_cb_type pu -net2as_overlap_frac 2e-1)
EOF

echo exit: $?
yq -io json '.Stdout |= from_yaml' $LOG
