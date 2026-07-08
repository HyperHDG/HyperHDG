#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
GRID=$DOMAIN-grid.geo.h5
MIKADO=$DOMAIN-mikado.geo.h5
LOG=$OUT/log.json
IMG=$OUT/$NAME.png
N=$((2**9+1))
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-10"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/ne18-14-fig8.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# h refers here to the edge length
ms="1 2 4 8 16 32"
for m in $ms; do
  tol=$(echo "scale=10; h=1/$N; ($m+.5)*h" | bc)
  python experiments/make_geo2.py --grid $N -o $OUT/domain-$m.geo.h5 \
    --dirichlet xmin=1 xmax=1 ymin=1 ymax=1 --dirichlet-tol $tol > /dev/null
done

parallel --results $LOG --progress --bar \
  "echo m: {1}; $BUILD/network -test diffusion -domain $OUT/domain-{1}.geo.h5 $KSP $NET -net2as_p {2}" \
  ::: $ms ::: 3 7 15 31
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-15-grid-dirichlet-tol-01.sh $LOG
