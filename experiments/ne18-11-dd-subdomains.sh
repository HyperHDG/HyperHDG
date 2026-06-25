#!/bin/bash
# Sweep the net2as "pu" number of subdomains (p*p) at a fixed overlap fraction. For every config we
# plot the residual history (-ksp_monitor_yaml) against wall-clock time. The monitor's time field
# already includes the factorization (it is logged at iteration 0), so the cost of more (smaller)
# subdomains is compared fairly against the iteration count it buys.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
: ${NP:=8}                  # MPI ranks (mat_cache is np-specific, keep fixed in a sweep)
: ${FRAC:=0.3}              # overlap distance / per-subdomain weighted diameter
: ${PS:="2 4 8 16"}         # subdomains per axis -> p*p subdomains
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
  echo "$BUILD/network -mat_cache $MAT -test constant -domain $DOMAIN $NET \
         -net2as_print_local -ksp_monitor_yaml -mem_max $@"
}

Q1="-net2as_cb_type q1 -net2as_cb_trim"
PU="-net2as_cb_type pu -net2as_overlap_frac .1"
for p in 1 2 4 8 16 32 64; do
  cmd $Q1 -net2as_p $p
  cmd $PU -net2as_p $p
done | parallel -j 1 --progress --bar --results $LOG

yq -io json '.Stdout |= from_yaml' $LOG
jq '.Stdout | {mem_max, t_ksp, t_iteration, iterations, label, p: .net2as.sz}' $LOG
