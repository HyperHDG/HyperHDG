#!/bin/bash
# Parallel-scaling / efficiency check: one fixed net2as pu config run at several MPI rank counts.
# Iterations need not be identical across rank counts (the graph partitioner is itself parallel, so
# the p*p subdomains can differ with np), but should be close; wall clock reveals parallel efficiency
# (subdomain factorizations are distributed across ranks but serial within a rank, and the coarse
# solve is a serial-ish bottleneck). Each rank count needs its own matrix cache (DoF numbering is
# np-specific). The same residual history (vs wall clock incl. factorization) is plotted per np.
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME/phd/nextcloud/networks/morgan-2026-05-20/net1/sca/}
: ${P:=4}                   # subdomains per axis -> p*p subdomains (fixed across rank counts)
: ${FRAC:=0.3}             # overlap distance / per-subdomain weighted diameter
: ${NPS:="1 2 4 8"}        # MPI rank counts to compare
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
LOG=$OUT/data.jsonl
PLOT=$OUT/$NAME.png
CUT=.4
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

# run <np>; per-np matrix cache (mat_cache is np-specific). archives raw json, appends residual points.
run() {
  np="$1"
  mpirun -n $np $BUILD/network -mat_cache $OUT/mat.$np.bin -test constant -domain $DOMAIN $NET \
         -net2as_cb_type pu -net2as_p $P -net2as_overlap_frac $FRAC -net2as_print_local -ksp_monitor_yaml 2>/dev/null \
    | yq -o json -I0 > $OUT/np$np.json
  jq -c --arg cfg "np=$np" '
      .iterations as $its |
      .ksp_monitor[] | {config: $cfg, t: .time, rnorm: .rnorm, iters: $its}' $OUT/np$np.json >> $LOG
}

: > $LOG
for np in $NPS; do run $np; done

echo "np  | iters | setup(it0)[s] | total[s] | sum_local_fac_nz | load_bal_mes"
for np in $NPS; do
  printf "%-3s | %5s | %13s | %8s | %16s | %s\n" "$np" \
    "$(jq -r '.iterations' $OUT/np$np.json)" \
    "$(jq -r '.ksp_monitor[0].time' $OUT/np$np.json)" \
    "$(jq -r '.ksp_monitor[-1].time' $OUT/np$np.json)" \
    "$(jq -r '[.net2as.local[].nz_fac]|add' $OUT/np$np.json)" \
    "$(jq -r '.net2as.load_bal.mes' $OUT/np$np.json)"
done

python experiments/plot.py -x t -y rnorm -g config --log y --marker "" \
       --xlabel "wall clock [s] (incl. factorization)" --ylabel "residual" \
       --title "$NAME: parallel scaling (pu p=$P, frac=$FRAC)" \
       --comment "$(git rev-parse --short HEAD)" --nshow --save $PLOT < $LOG
sxiv $PLOT
