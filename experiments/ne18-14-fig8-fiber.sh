#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
LOG=$OUT/log.json
IMG=$OUT/$NAME.png
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET
cp $BUILD/network experiments/make_geo2.py experiments/ne18-14-fig8-plot.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# gortz.pdf fig. 8 (center panel) / table 2 "Fiber": heat conductivity on the uniform
# random fiber network of sec. 6.1 (mass 1000, fiber length r = 0.05, ~0.36M nodes),
# constant weight gamma = 1, f = M1, u = 0 on the boundary nodes (the clipped fiber
# endpoints on the unit-square edge; dirichlet-tol tight so only those are pinned).
# Paper reference rates (tau_avg, tau_max): H=1/4: (0.29,0.48), 1/8: (0.33,0.43),
# 1/16: (0.39,0.49), 1/32: (0.42,0.47).
python experiments/make_geo2.py --mikado 1000 -o $DOMAIN \
  --dirichlet xmin=1 xmax=1 ymin=1 ymax=1 --dirichlet-tol 1e-5

# net2as coarse grid has p+1 cells per axis: p = H^-1 - 1 gives the paper's H = 1/4..1/32
for p in 3 7 15 31; do
  { echo "Hinv: $((p+1))";
    mpirun -n 1 $BUILD/network -test diffusion -domain $DOMAIN $NET -net2as_p $p \
      -ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-10 -ksp_norm_type unpreconditioned; } | yq -o json -I0 >> $LOG
done

experiments/ne18-14-fig8-plot.sh $LOG $IMG
