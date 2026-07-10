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
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
IMG=$OUT/$NAME.png
# conditioning of ONE implicit timowave step vs timestep length, on the ne18-23 grid
# (unit square 129, unit props; fundamental period T1 ~ 2.1 measured by ne18-23 --
# the dispersion estimate 1.45 misses the sqrt(2) grid-orientation factor: both edge
# families carry mass, only the aligned one carries shear).  The theta scheme's local
# operator carries the mass shift C_u/(theta dt^2): halving dt quadruples the shift,
# so smaller steps push more of the spectrum into the mass-dominated (DD-friendly)
# regime.  The shift passes the fundamental lambda_1 = omega_1^2 ~ 9 around
# dt ~ T1/3: most of the ladder T1/2 .. T1/32 sits beyond it, progressively burying
# the low (DD-hard) end of the spectrum under the mass term.
# Metric: energy-norm error vs CG iteration (the fig8-style curves), one curve per dt,
# at p = 15 i.e. H = 1/16.  enorm's reference solve happens inside the first timestep
# (timowave calls KSPMonitorYAML_Setup there; meaningful because -nt 1).
: ${THETA:=0.5}
: ${DEG:=3}
: ${P:=15}
T1=2.1
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

DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 1e-3"
$PYTHON experiments/make_geo2.py --grid 129 $DIR -o $GRID | tee $LOGG

for f in 2 4 8 16 32; do
  DT=$($PYTHON -c "print($T1/$f)")
  { echo "dtfrac: T1/$f";
    $MPIRUN -n $NP $BUILD/timowave -test constant -domain $GRID -deg $DEG \
      -theta $THETA -nt 1 -T $DT $KSP $NET; } | yq -o json -I0 >> $LOG
done

# summary: iterations of the single implicit solve per dt
jq -r '[.dtfrac, .dt, .iterations, .ref_its] | @tsv' $LOG | column -t

# fig8-style energy-error curves, one per timestep length
jq -c '{dtfrac, dt} + (.ksp_monitor[] | {it, enorm})' $LOG \
  | $PYTHON experiments/plot.py -x it -y enorm -g dtfrac --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG --tikz ${IMG%.png}
