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

# run <label> <extra net2as args...>; archives raw json in $OUT and appends one json line per
# residual-history point. t = ksp_monitor wall clock (already includes the factorization at it=0).
run() {
  label="$1"; shift
  mpirun -n $NP $BUILD/network -mat_cache $MAT -test constant -domain $DOMAIN $NET \
         -net2as_cb_type pu -net2as_overlap_frac $FRAC -net2as_print_local -ksp_monitor_yaml "$@" 2>/dev/null \
    | yq -o json -I0 > "$OUT/$label.json"
  jq -c --arg cfg "$label" '
      .iterations as $its |
      .ksp_monitor[] | {config: $cfg, t: .time, rnorm: .rnorm, iters: $its}' "$OUT/$label.json" >> $LOG
}

: > $LOG
for p in $PS; do run "p$p" -net2as_p $p; done

echo "config | nsub | iters | wall_end[s] | sum_fac_nz | max_fac_nz | sum_fac_t[s]"
for p in $PS; do
  printf "p=%-3s | %4s | %5s | %11s | %10s | %10s | %s\n" "$p" "$((p*p))" \
    "$(jq -r '.iterations' "$OUT/p$p.json")" \
    "$(jq -r '.ksp_monitor[-1].time' "$OUT/p$p.json")" \
    "$(jq -r '[.net2as.local[].nz_fac]|add' "$OUT/p$p.json")" \
    "$(jq -r '[.net2as.local[].nz_fac]|max' "$OUT/p$p.json")" \
    "$(jq -r '[.net2as.local[].time]|add' "$OUT/p$p.json")"
done

python experiments/plot.py -x t -y rnorm -g config --log y --marker "" \
       --xlabel "wall clock [s] (incl. factorization)" --ylabel "residual" \
       --title "$NAME: convergence vs #subdomains (overlap frac=$FRAC)" \
       --comment "$(git rev-parse --short HEAD)" --nshow --save $PLOT < $LOG
sxiv $PLOT
