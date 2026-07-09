#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
FIBER2RAW=$DOMAIN-fiber2raw.geo.h5
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

# ne18-20 with the two lessons of the solution inspection folded in:
#  - RAW properties, no --prop-cutoff: with the honest RHS (welds unloaded, see below) the
#    net2as iteration counts are INDEPENDENT of the cutoff level (quarter probe: c5/c2/
#    c0.5/raw identical, 23/44/55 its at p=2/7/13) - the exact local solves absorb the
#    full 1e10..1e20 weld contrast, so the cutoff regularisation is unnecessary.
#  - the 'constant' test now opts into massless_unloaded: virtual weld edges (mass 0,
#    fiber_id -1, 34% of edges) carry no body load, which removes the huge localized
#    weld-chain displacements (max|u| 32 -> 0.77 on the quarter) and makes the solution
#    a clean z-dominated sag.
# Rotation rigidities G_xI_x, E_1I_1, E_2I_2 x 1e6 keep l_c ~ 1700 um >= R_max = 1000 um
# (shear-dominated regime, cf. the gortz_constants l_c section).  p = 3..31 = fig8's
# H^-1 = 4..32.  ne18-20's p=31 run hit DIVERGED_NANORINF with floored props; watch
# whether raw props reproduce it.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2RAW $DIR $SUB $REG | tee $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2RAW --cells 4 8 16 32 | tee $OUT/gortz-fiber2raw.txt

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {3} -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {2}; echo domain: {1}" \
  ::: fiber2raw ::: 3 7 15 31 ::: constant
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
