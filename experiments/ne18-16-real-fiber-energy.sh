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
FIBER1=$DOMAIN-fiber1.geo.h5
FIBER2=$DOMAIN-fiber2.geo.h5
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
IMG=$OUT/$NAME.png
N=513
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-9"
NP=$(nproc)
MPIRUN=spack/$PRESET/.spack-env/view/bin/mpirun
GORTZ=$OUT/gortz.txt

export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/ne18-14-fig8-plot.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
PROP="--prop-cutoff 5"

#python experiments/make_geo2.py --grid $((2**9+1)) -o $GRID $DIR --dirichlet-tol 1e-3 >> $LOGG
#python experiments/make_geo2.py --mikado 1000 -o $MIKADO $DIR >> $LOGG
#python experiments/make_geo2.py -i domains/fiber-2026-05-20/net1/sca -o $FIBER1 $DIR $SUB --no-props | tee $LOGG
python experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2 $DIR $SUB --no-props | tee $LOGG

python experiments/gortz_constants.py $FIBER2 --mu --cells 3 4 8 10 12 14 | tee $GORTZ

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {3} -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {2}; echo domain: {1}" \
  ::: fiber2 ::: 2 3 7 11 13 ::: constant
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
