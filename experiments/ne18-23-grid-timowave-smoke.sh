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
WAVE=$OUT/wave.vtkhdf
# collective vtkhdf write deadlocks on NFS locking (see ne18-22) - node-local scratch
WAVE_SCRATCH=$(mktemp -d ${TMPDIR:-/tmp}/$NAME.XXXXXX)/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png
# ne18-22's fiber smoke on the clean homogeneous baseline: unit grid, unit props
# (l_c = 1, c_s = 1, mass = edge length -> mu = 1).  Timoshenko dispersion at
# k = pi*sqrt2/extent gives omega_1 = 4.334 -> T1 = 1.45, but on the grid only the
# aligned edge family carries shear for a given gradient while BOTH carry mass:
# c_eff = c_s/sqrt(2), omega_1 = 3.14 -> T1 = 2.05; measured (first run of this
# script, sag turnaround at t ~ 1.05) T1 ~ 2.1.  One measured period.
: ${NT:=20}
: ${T:=2.1}
: ${THETA:=0.5}
: ${DEG:=3}
NET="-pc_type net2as -net2as_p 7 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
KSP="-ksp_monitor_yaml -ksp_rtol 1e-9"
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
   experiments/netvis.py experiments/energy.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-19's Timoshenko grid convention: unit square, full clamp (type 63) on all four
# borders, tol below h = 1/128 so only the boundary ring is pinned.  Default (unit)
# properties kept so gortz_constants can print the wave speeds / T1 cross-check;
# mass = edge length means density 1, same coefficients the solver defaults to.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 1e-3"

$PYTHON experiments/make_geo2.py --grid 129 $DIR -o $GRID | tee $LOGG

$PYTHON experiments/gortz_constants.py $GRID --cells 4 8 16 \
  --csv $OUT/gortz-grid.csv | tee $OUT/gortz-grid.txt

stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test constant -domain $GRID -deg $DEG \
  -theta $THETA -nt $NT -T $T $KSP $NET -plot $WAVE_SCRATCH -print_timestep | tee $LOG
mv $WAVE_SCRATCH $WAVE
rmdir $(dirname $WAVE_SCRATCH)

# energy exchange over time: kinetic <-> strain shows the oscillation quantitatively
$PYTHON experiments/energy.py $WAVE -o $NRG || echo "energy plot failed (non-fatal)"

# renders (need pvpython): static sag ~0.14 on the unit square, dynamic max ~2x -> warp 1
# is already a ~30% deflection.  Still: 4 overlaid steps ~ 0, T1/4, T1/2 (max sag), 3T1/4.
experiments/netvis.py $WAVE --beams 1 --view iso --warp-scale 1 --show 0 --duration 8 -o $VID \
  || echo "netvis animation failed (non-fatal, e.g. no pvpython)"
experiments/netvis.py $WAVE --beams 1 --view pside --warp-scale 1 --show 0 --axis 0 \
  --frames 0,0.5,1.05,1.6 --frame-colors viridis -o $IMG \
  || echo "netvis still failed (non-fatal, e.g. no pvpython)"
