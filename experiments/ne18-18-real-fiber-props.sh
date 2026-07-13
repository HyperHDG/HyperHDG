#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
FIBER2P=$DOMAIN-fiber2p.geo.h5
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

# fiber2 with TRUE properties (kg/um/s; ribbons ~21x1.5 um) instead of ne18-16's
# --no-props.  --prop-cutoff 5 floors the floppy tail (EA/GA contrast 1e10 -> 4e2,
# EI/GI 1e20 -> 3e5); cutoff 2 vs 5 measured identical, so contrast is not the binding
# constraint.  The binding constraint is the bending length l_c = sqrt(EI/GA) ~ 1.7 um
# (thin direction) << R = 285..1333 um: the shear term GA|u' - r x t|^2 penalizes the
# componentwise q1 interpolant by (R/l_c)^2 relative to the bending modes it must
# approximate, so the coarse level is inert and the method is one-level (quarter-crop
# probe: its = 8n, kappa = 1.6 n^2).  Arm fiber2reg scales the rotation-block
# stiffnesses G_xI_x, E_1I_1, E_2I_2 by s = 1e6 (l_c *= 1e3, ties r to u), putting all
# R <= ext/(2*3) = 1333 <= l_c_p50 = 1700: probe shows saturation already at s where
# l_c ~ R_max/2, so s = 1e6 is safely in the flat regime.  gortz_constants now prints
# the l_c section for both files.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
CUT="--prop-cutoff 5"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2P $DIR $SUB $CUT | tee $LOGG
$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2REG $DIR $SUB $CUT $REG | tee -a $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2P --cells 3 4 8 12 14 | tee $OUT/gortz-fiber2p.txt
$PYTHON experiments/gortz_constants.py $FIBER2REG --cells 3 4 8 12 14 | tee $OUT/gortz-fiber2reg.txt

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {3} -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {2}; echo domain: {1}" \
  ::: fiber2p fiber2reg ::: 2 3 7 11 13 ::: constant
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
