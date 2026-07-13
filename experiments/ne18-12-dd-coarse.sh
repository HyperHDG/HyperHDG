#!/bin/bash
# Compare net2as pu coarse-space richness at a fixed overlap: flat PU vs linear {1,x,y} vs bilinear
# {1,x,y,xy} (q1-like, not full quadratic), against the q1-trim baseline (same p*p subdomains). Tests
# whether richer coarse enrichment reduces pu iterations. Plots residual vs wall clock (which already
# includes factorization at iteration 0); each config's raw json is archived under $OUT.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
: ${NP:=8}                  # MPI ranks (mat_cache is np-specific, keep fixed in a sweep)
: ${P:=4}                   # subdomains per axis -> p*p subdomains
: ${FRAC:=0.05}             # overlap distance / graph diameter (small, where the coarse space matters most)
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
MAT=$OUT/mat.bin
LOG=$OUT/data.jsonl
PLOT=$OUT/$NAME.png
CUT=.4
# unpreconditioned norm so every config converges in the same true residual ||b-Ax||
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
    | yq -o json -I0 > $OUT/$label.json
  jq -c --arg cfg "$label" '
      .iterations as $its |
      .ksp_monitor[] | {config: $cfg, t: .time, rnorm: .rnorm, iters: $its}' $OUT/$label.json >> $LOG
}

: > $LOG
run "q1trim"   -net2as_cb_type q1 -net2as_cb_trim
run "flat"     -net2as_cb_type pu -net2as_overlap_frac $FRAC
run "linear"   -net2as_cb_type pu -net2as_overlap_frac $FRAC -net2as_pux_dim 2
run "bilinear" -net2as_cb_type pu -net2as_overlap_frac $FRAC -net2as_pux_dim 2 -net2as_pu_xy

echo "config    | iters | wall_end[s] | coarse_size"
for f in q1trim flat linear bilinear; do
  printf "%-9s | %5s | %11s | %s\n" "$f" \
    "$(jq -r '.iterations' $OUT/$f.json)" \
    "$(jq -r '.ksp_monitor[-1].time' $OUT/$f.json)" \
    "$(jq -r '.net2as.coarse.size' $OUT/$f.json)"
done

python experiments/plot.py -x t -y rnorm -g config --log y --marker "" \
       --xlabel "wall clock [s] (incl. factorization)" --ylabel "residual" \
       --title "$NAME: coarse-space richness (p=$P, frac=$FRAC)" \
       --comment "$(git rev-parse --short HEAD)" --nshow --save $PLOT < $LOG
sxiv $PLOT
