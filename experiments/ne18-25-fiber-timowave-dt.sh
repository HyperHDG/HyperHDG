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
# ne18-24's single-step dt ladder moved to the quarter fiber probe (ne18-22 domain:
# fiber2 raw props, welds unloaded via the constant test, rotation rigidities x1e6,
# --clamp-xy .25).  On the homogeneous grid the ladder was FLAT (16-17 its, q1
# two-level already optimal).  Here the preconditioner's weakness IS the global
# low-energy modes (kappa grows with H^-1, ne18-21), which are exactly what the mass
# shift C_u/(theta dt^2) suppresses -- so the dt-dependence should show, if anywhere.
# T1 = 2.0e-6 s is the sag-mode period measured by ne18-22 at deg 3 (dispersion
# estimate 1.7e-6).  Metric: energy-norm error vs CG iteration, one curve per dt,
# p = 15 i.e. H = 1/16.  enorm's reference solve runs inside the first timestep
# (meaningful because -nt 1); note the trace matrix depends on dt through the local
# condensation, so nothing is reusable across the ladder.
: ${THETA:=0.5}
: ${DEG:=3}
: ${P:=15}
T1=2.0e-6
NET="-pc_type net2as -net2as_p $P -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
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

for f in 2 4 8 16 32; do
  DT=$($PYTHON -c "print($T1/$f)")
  { echo "dtfrac: T1/$f";
    $MPIRUN -n $NP $BUILD/timowave -test constant -domain $FIBER2RAWQ -deg $DEG \
      -theta $THETA -nt 1 -T $DT $KSP $NET; } | yq -o json -I0 >> $LOG
done

# summary: iterations of the single implicit solve per dt
jq -r '[.dtfrac, .dt, .iterations, .ref_its] | @tsv' $LOG | column -t

# fig8-style energy-error curves, one per timestep length
jq -c '{dtfrac, dt} + (.ksp_monitor[] | {it, enorm})' $LOG \
  | $PYTHON experiments/plot.py -x it -y enorm -g dtfrac --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG --tikz ${IMG%.png}
