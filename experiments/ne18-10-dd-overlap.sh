#!/bin/bash
# Sweep the net2as "pu" overlap distance (as a fraction of the weighted graph diameter) at a fixed
# number of subdomains, plus the q1-trim baseline (same subdomain count). For every config we plot
# the residual history (-ksp_monitor_yaml) against wall time, where wall time = factorization time
# + iteration time, so configurations that trade a cheaper factorization for more iterations are
# compared fairly. q1-trim is included because it produces p*p subdomains, matching pu.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
: ${NP:=8}                              # MPI ranks (mat_cache is np-specific, keep fixed in a sweep)
: ${P:=4}                               # subdomains per axis -> p*p subdomains
: ${FRACS:="0.04 0.05 0.075 0.10 0.15"} # overlap distance / graph diameter
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
MAT=$OUT/mat.bin
LOG=$OUT/data.jsonl
PLOT=$OUT/$NAME.png
CUT=.4
NET="-pc_type net2as -net2as_pc_factor_mat_solver_type mumps"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/{make_geo2,plot}.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py -i $INPUT --clamp-xy $CUT -o $DOMAIN \
       --dirichlet xmax=63 xmin=63 ymin=63 ymax=63

# run <config-label> <extra net2as args...>; appends one json line per residual-history point.
# t = (mean per-rank local factorization + coarse factorization) + (iteration wall time, rebased to 0)
run() {
  label="$1"; shift
  mpirun -n $NP $BUILD/network -mat_cache $MAT -test constant -domain $DOMAIN $NET \
         -net2as_p $P -net2as_print_local -ksp_monitor_yaml "$@" 2>/dev/null \
    | yq -o json -I0 \
    | jq -c --arg cfg "$label" '
        (([.net2as.local[].time] | add) / .mpi.sz + (.net2as.coarse.time // 0)) as $tfac |
        .ksp_monitor[0].time as $t0 |
        .iterations as $its |
        .ksp_monitor[] | {config: $cfg, t: ($tfac + (.time - $t0)), rnorm: .rnorm, iters: $its}'
}

: > $LOG
run "q1trim" -net2as_cb_type q1 -net2as_cb_trim >> $LOG
for f in $FRACS; do
  run "f=$f" -net2as_cb_type pu -net2as_overlap_frac $f >> $LOG
done

echo "config | iters | wall_end[s] (fac+iter)"
jq -rs 'group_by(.config)[] | "\(.[0].config) | \(.[0].iters) | \(([.[].t] | max * 1000 | round / 1000))"' $LOG

python experiments/plot.py -x t -y rnorm -g config --log y --marker "" \
       --xlabel "wall time: factorization + iterations [s]" --ylabel "residual" \
       --title "$NAME: convergence vs overlap (p=$P, $((P*P)) subdomains)" \
       --comment "$(git rev-parse --short HEAD)" --nshow --save $PLOT < $LOG
sxiv $PLOT
