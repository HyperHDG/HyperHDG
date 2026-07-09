#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
FIBER2REG=$DOMAIN-fiber2reg.geo.h5
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
IMG=$OUT/$NAME.png
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-9"
NP=$(nproc)
MPIRUN=spack/$PRESET/.spack-env/view/bin/mpirun
[ -x "$MPIRUN" ] || MPIRUN=mpirun
# spack env python when usable (working matplotlib on the pde cluster), else system python
if [ -z "${PYTHON:-}" ]; then
  PYTHON=spack/$PRESET/.spack-env/view/bin/python
  $PYTHON -c 'import numpy, h5py, pandas, matplotlib' >/dev/null 2>&1 || PYTHON=python
fi
export PYTHON

export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/ne18-14-fig8-01-residual-iteration.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-18's regularised arm (fiber2 with true props, floppy tail floored at the 5th
# percentile, rotation rigidities G_xI_x, E_1I_1, E_2I_2 x 1e6 so l_c ~ 1700 um >= every
# subdomain radius) extended to the fig8 subdomain range p = 3..31, i.e. H^-1 = 4..32
# (H = 2000..250 um).  Note l_c_p50 = 1700 sits between R(p=3) = 1000 and R(p=7) = 500,
# so the whole sweep stays in the shear-dominated regime.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
CUT="--prop-cutoff 5"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2REG $DIR $SUB $CUT $REG | tee $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2REG --cells 4 8 16 32 | tee $OUT/gortz-fiber2reg.txt

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {3} -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {2}; echo domain: {1}" \
  ::: fiber2reg ::: 3 7 15 31 ::: constant
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
