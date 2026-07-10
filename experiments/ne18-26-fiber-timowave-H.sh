#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
FIBER2RAWQ=$DOMAIN-fiber2rawq.geo.h5
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
IMG=$OUT/$NAME.png
# the H-sweep companion to ne18-25's dt ladder: ONE implicit timowave step at the
# longest physical timestep dt = T1/2 (T1 = 2.0e-6, the sag period measured by
# ne18-22 at deg 3), subdomain count p = 3 7 15 31 i.e. H^-1 = 4 8 16 32 -- the wave
# analogue of the ne18-21 stationary sweep on the quarter probe.  At dt = T1/2 the
# mass shift is minimal (ne18-25: it buys only ~-14% by T1/32), so expectation is
# stationary-like its growth with H^-1; the p=15 rung must reproduce ne18-25's
# T1/2 run (64 its) exactly.
: ${THETA:=0.5}
: ${DEG:=3}
T1=2.0e-6
DT=$(python -c "print($T1/2)" 2>/dev/null || echo 1.0e-6)
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-9"
: ${NP:=$(nproc)}
MPIRUN=spack/$PRESET/.spack-env/view/bin/mpirun
[ -x "$MPIRUN" ] || MPIRUN=mpirun
# spack env python when usable (working matplotlib on the pde cluster), else system python
if [ -z "${PYTHON:-}" ]; then
  PYTHON=spack/$PRESET/.spack-env/view/bin/python
  $PYTHON -c 'import numpy, h5py, pandas, matplotlib' >/dev/null 2>&1 || PYTHON=python
fi
export PYTHON

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-22's domain build (= ne18-21 cut to the quarter probe)
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"
CLAMP="--clamp-xy .25"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2RAWQ \
  $DIR $SUB $REG $CLAMP | tee $LOGG

for p in 3 7 15 31; do
  { echo "H: 1/$((p+1))";
    $MPIRUN -n $NP $BUILD/timowave -test constant -domain $FIBER2RAWQ -deg $DEG \
      -theta $THETA -nt 1 -T $DT $KSP $NET -net2as_p $p; } | yq -o json -I0 >> $LOG
done

# summary: iterations of the single implicit solve per coarse scale
jq -r '[.H, .iterations, .ref_its] | @tsv' $LOG | column -t

# fig8-style energy-error curves, one per H
jq -c '{H, dt} + (.ksp_monitor[] | {it, enorm})' $LOG \
  | $PYTHON experiments/plot.py -x it -y enorm -g H --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG
