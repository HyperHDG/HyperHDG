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

# gortz.pdf fig. 8 (left panel) / table 2 "Grid": heat conductivity on a uniform regular
# grid with (2^9+1)^2 nodes, constant coefficient 1, f = M1 (constant source per unit
# length), u = 0 on the whole boundary; subspace decomposition preconditioner with coarse
# scale H = 1/p for p subdomains per axis.
# dirichlet-tol < h = 1/512 pins ONLY the geometric boundary ring (the paper's u = 0 on
# the boundary nodes); the default 2e-2 would pin ~10 node layers, moving the effective
# boundary off the coarse grid lines and degrading the near-boundary decomposition
python experiments/make_geo2.py --grid 513 -o $DOMAIN \
  --dirichlet xmin=1 xmax=1 ymin=1 ymax=1 --dirichlet-tol 1e-3

# net2as builds the coarse Q1 grid with p+1 cells per axis, so p = H^-1 - 1 gives the paper's
# H = 1/4..1/32 with the coarse mesh aligned to the 512-cell network grid (H multiple of h);
# cb_trim = BC-conforming coarse space Q_H, subdomain cover keeps all patches
for p in 3 7 15 31; do
  { echo "Hinv: $((p+1))";
    mpirun -n 1 $BUILD/network -test diffusion -domain $DOMAIN $NET -net2as_p $p \
      -ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-10 -ksp_norm_type unpreconditioned; } | yq -o json -I0 >> $LOG
done

experiments/ne18-14-fig8-plot.sh $LOG $IMG
