#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOM=$OUT/domain.geo.h5
LOG=$OUT/log.json
IMG=$OUT/$NAME.png
N=513
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-10"
NP=$(nproc)
MPIRUN=spack/$PRESET/.spack-env/view/bin/mpirun
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/ne18-14-fig8*.sh experiments/plot.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

python experiments/make_geo2.py -i 'domains/fiber-2026-05-20/net1/sca' -o $DOM \
  --dirichlet xmin=1 xmax=1 ymin=1 ymax=1 --dirichlet-tol 2e-2

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test diffusion -domain {1} $KSP $NET -net2as_p {2}" \
  ::: $DOM ::: 3 7 15 31
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
