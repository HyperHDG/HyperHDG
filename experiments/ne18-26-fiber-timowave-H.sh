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
# the H-sweep companion to ne18-25's dt ladder, on the FULL fiber2 domain (ne18-21's
# canonical config) so every rung H = 8000/(p+1) = 2000..250 um stays well above
# R0 = 62.5 um and the Goertz constants do not degrade (a quarter-domain attempt put
# H(p=31) = R0 exactly and its blew up to 112 -- scale effect, not the wave operator).
# ONE implicit timowave step at the longest physical timestep dt = T1/2.  T1 = 8.0e-6 s:
# 4x the sag period measured on the quarter (ne18-22, shear regime T ~ L/c_s); the
# dispersion estimate at the full-domain k has l_c*k ~ 0.95, i.e. the x1e6
# regularisation puts the crossover right at the domain scale (T1 ~ 9.6e-6) -- the
# factor is immaterial for its (ne18-25: -14% over a factor 16 in dt).
# Second arm: the STATIONARY sweep (network, same test/domain/rungs) = ne18-21
# reproduced in-run, for an apples-to-apples wave-vs-stationary comparison.
: ${THETA:=0.5}
: ${DEG:=3}
T1=8.0e-6
DT=$(python -c "print($T1/2)" 2>/dev/null || echo 4.0e-6)
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
cmake --build --preset $PRESET --target timowave --target network
cp $BUILD/timowave $BUILD/network experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-21's full-domain build
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2RAW \
  $DIR $SUB $REG | tee $LOGG

for p in 3 7 15 31; do
  { echo "H: 1/$((p+1))"; echo "solver: wave";
    $MPIRUN -n $NP $BUILD/timowave -test constant -domain $FIBER2RAW -deg $DEG \
      -theta $THETA -nt 1 -T $DT $KSP $NET -net2as_p $p; } | yq -o json -I0 >> $LOG
  { echo "H: 1/$((p+1))"; echo "solver: stat";
    $MPIRUN -n $NP $BUILD/network -test constant -domain $FIBER2RAW \
      $KSP $NET -net2as_p $p; } | yq -o json -I0 >> $LOG
done

# summary: iterations of the single implicit solve per coarse scale and arm
jq -r '[.H, .solver, .iterations, .ref_its] | @tsv' $LOG | column -t

# fig8-style energy-error curves, one per H; one figure per arm (wave / stat)
jq -c '{H, solver} + (.ksp_monitor[] | {it, enorm})' $LOG \
  | $PYTHON experiments/plot.py -x it -y enorm -g H --group0 solver --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG --tikz ${IMG%.png}
