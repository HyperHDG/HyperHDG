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
: ${NP:=1}                              # MPI ranks (mat_cache is np-specific, keep fixed in a sweep)
: ${P:=4}                               # subdomains per axis -> p*p subdomains
: ${FRACS:="0.1 0.2 0.3"}      # overlap distance / per-subdomain weighted diameter
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
MAT=$OUT/mat.bin
LOG=$OUT/data.jsonl
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

# run <config-label> <extra net2as args...>; archives raw json in $OUT and appends one json line per
# residual-history point. t = ksp_monitor wall clock (already includes the factorization at it=0).
run() {
  label="$1"; shift
  mpirun -n $NP $BUILD/network -mat_cache $MAT -test constant -domain $DOMAIN $NET \
         -net2as_p $P -net2as_print_local -ksp_monitor_yaml "$@" 2>/dev/null \
    | yq -o json -I0 > "$OUT/$label.json"
  jq -c --arg cfg "$label" '
      .iterations as $its |
      .ksp_monitor[] | {config: $cfg, t: .time, rnorm: .rnorm, iters: $its}' "$OUT/$label.json" >> $LOG
}

: > $LOG
LABELS="q1trim"
run "q1trim" -net2as_cb_type q1 -net2as_cb_trim
for f in $FRACS; do
  run "f=$f" -net2as_cb_type pu -net2as_overlap_frac $f
  LABELS="$LABELS f=$f"
done

# total/max subdomain factorization nnz and summed factor time show the small-subdomain advantage
echo "config   | iters | wall_end[s] | sum_fac_nz | max_fac_nz | sum_fac_t[s]"
for l in $LABELS; do
  printf "%-8s | %5s | %11s | %10s | %10s | %s\n" "$l" \
    "$(jq -r '.iterations' "$OUT/$l.json")" \
    "$(jq -r '.ksp_monitor[-1].time' "$OUT/$l.json")" \
    "$(jq -r '[.net2as.local[].nz_fac]|add' "$OUT/$l.json")" \
    "$(jq -r '[.net2as.local[].nz_fac]|max' "$OUT/$l.json")" \
    "$(jq -r '[.net2as.local[].time]|add' "$OUT/$l.json")"
done

python experiments/plot.py -x t -y rnorm -g config --log y --marker "" \
       --xlabel "wall clock [s] (incl. factorization)" --ylabel "residual" \
       --title "$NAME: convergence vs overlap (p=$P, $((P*P)) subdomains)" \
       --comment "$(git rev-parse --short HEAD)" --nshow --save $PLOT < $LOG
sxiv $PLOT
