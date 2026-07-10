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
WAVE=$OUT/wave.vtkhdf
# the collective (MPI-IO) vtkhdf write deadlocks on NFS byte-range locking (ranks stuck
# in nlmclnt_wait on pde12's NFS home) - write to node-local scratch, then move to $OUT.
# Fine on a single node; multi-node runs need a real parallel FS here instead.
WAVE_SCRATCH=$(mktemp -d ${TMPDIR:-/tmp}/$NAME.XXXXXX)/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png
# time horizon: the quarter probe has xy extent ~2000 um and transverse (shear-regime)
# speed c_s = sqrt(kGA/(m/l)) ~ 1.7e9 um/s (gortz_constants wave-speed section), so the
# fundamental period is T1 ~ 2*extent/c_s ~ 2.3e-6 s -- microseconds, not milliseconds
# (measured at deg 1: global z-sag half-period ~4.7e-6, but deg 1 is ~11x too soft, see
# below).  Default: ~2 periods in 20 steps, enough to see the network sag under the
# constant force and swing back.
# theta: weight on the new time level of the one-step theta scheme; stiff modes amplify
# by -(1-theta)/theta per step, so theta < 0.5 blows up (theta 0.25 grew x3/step here,
# on the tiny grid, everywhere -- the old "0.25 = CN" help text was wrong).  0.5 = CN.
# deg: the deg-1 static limit came out ~11x softer than the network (deg 3) static
# solve (max|u| 8.90 vs 0.754 um, relaxation vs elliptic); hybrid energy >> physical
# says trace and bulk barely agree at deg 1.  Run the smoke at deg 3.
: ${NT:=20}
: ${T:=5e-6}
: ${THETA:=0.5}
: ${DEG:=3}
# coarse LU instead of the default Cholesky, same reason as ne18-21 (numerically
# indefinite coarse matrix under the x1e6 rotation rigidities at high p)
NET="-pc_type net2as -net2as_p 7 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
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
# HDF5's own flock on NFS (h5py readers, serial writes); the MPI-IO locking above is
# ROMIO's and needs the scratch redirect regardless
export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/netvis.py experiments/energy.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-21's canonical stationary config (fiber2, RAW props, welds unloaded via the
# constant test's massless_unloaded, rotation rigidities x1e6) cut down to the quarter
# probe for a first timowave smoke test: the constant body force switches on at t=0
# with the network at rest, so it sags towards the ne18-21 static solution (max|u|
# ~ 0.77 um on the quarter) and oscillates about it.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"
CLAMP="--clamp-xy .25"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2RAWQ \
  $DIR $SUB $REG $CLAMP | tee $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2RAWQ --cells 4 8 16 \
  --csv $OUT/gortz-fiber2rawq.csv | tee $OUT/gortz-fiber2rawq.txt

# stdbuf: keep PetscPrintf progress line-buffered through the tee pipe
stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test constant -domain $FIBER2RAWQ -deg $DEG \
  -theta $THETA -nt $NT -T $T $KSP $NET -plot $WAVE_SCRATCH -print_timestep | tee $LOG
mv $WAVE_SCRATCH $WAVE
rmdir $(dirname $WAVE_SCRATCH)

# energy exchange over time: kinetic <-> strain shows the oscillation quantitatively
$PYTHON experiments/energy.py $WAVE -o $NRG || echo "energy plot failed (non-fatal)"

# renders (need pvpython): max dynamic |u| ~ 1.5 um on the 2000 um probe, warp x100 for
# a ~8% visible deflection.  Still: 4 overlaid steps ~ 0, T1/2 (max sag), T1, 3T1/2.
experiments/netvis.py $WAVE --beams 1 --view iso --warp-scale 100 --show 0 -o $VID \
  || echo "netvis animation failed (non-fatal, e.g. no pvpython)"
experiments/netvis.py $WAVE --beams 1 --view pside --warp-scale 100 --show 0 --axis 0 \
  --frames 0,1.25e-6,2.5e-6,3.75e-6 --frame-colors viridis -o $IMG \
  || echo "netvis still failed (non-fatal, e.g. no pvpython)"
