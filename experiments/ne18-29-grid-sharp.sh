#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
GRID=$OUT/grid.geo.h5
WAVE=$OUT/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png

# SHARP-tap variant of ne18-27: same 500^2 grid, but std_x = h (tap sampled by ~2
# nodes) so the pulse pumps near-Nyquist lattice modes.  ne18-27 measured a perfectly
# isotropic front at c_eff = c_beam/sqrt(2) = 0.71 -- at wavelengths ~3h the grid is an
# effective isotropic plate (per-direction stiffness from one beam family, inertia from
# both).  Lattice anisotropy (axis vs diagonal split) can only appear near lattice-scale
# wavelengths; this run probes exactly that with energy_polar's r_q(theta).
# dt = T/NT = 2.5e-4 <= std_t/2.7 resolves the sharp pulse; T = 0.2 puts the front at
# r ~ 0.14 (~70 cells), enough for a few-percent speed split to separate.
: ${N:=500}
: ${NT:=800}
: ${T:=0.2}
: ${THETA:=0.5}
: ${DEG:=2}
: ${STD_X:=2e-3}
: ${STD_T:=6.7e-4}
: ${TAP:=3.3e-3}
: ${ENERGY:=5e-4}
# new plot granularity flags (ne18-27ff): stride 10 -> 61 stored frames; values disp ->
# 3 of 18 components (renderers must then address the filtered layout: --warp-by
# values:0,1,2 --color-by values:2); per-cell energies kept at this scale for energy.py
# and for validating the upcoming in-situ histograms against a full-fidelity file.
: ${STRIDE:=10}
: ${PLOT_VALUES:=disp}
: ${PLOT_PROPS:=none}

NET="-pc_type net2as -net2as_p 7 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
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
# collective vtkhdf writes onto NFS (cluster home2) deadlock in NLM byte-range locking;
# ompio's skip-locking is safe single-node and verified correct vs a local-disk reference
# (2026-07-11: hang without it, output equal to solver run-to-run noise ~3e-11 with it).
# Harmless on local disks, so exported unconditionally.
export OMPI_MCA_fs_ufs_lock_algorithm=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/netvis.py experiments/energy.py experiments/energy_rt.py experiments/energy_polar.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

$PYTHON experiments/make_geo2.py --grid $N \
  --dirichlet xmin=63 xmax=63 ymin=63 ymax=63 -o $GRID | tee $LOGG

$PYTHON experiments/gortz_constants.py $GRID --mu --cells 4 8 16 \
  --csv $OUT/gortz-grid.csv | tee $OUT/gortz-grid.txt

stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test drumhead -domain $GRID -deg $DEG \
  -theta $THETA -nt $NT -T $T \
  -std_x $STD_X -std_t $STD_T -tap_time $TAP -energy $ENERGY \
  -plot_stride $STRIDE -plot_values $PLOT_VALUES -plot_props $PLOT_PROPS \
  $KSP $NET -plot $WAVE -print_timestep | tee $LOG

$PYTHON experiments/energy.py $WAVE -o $NRG || echo "energy plot failed (non-fatal)"
$PYTHON experiments/energy_rt.py $WAVE -o $OUT/$NAME-rt.png \
  || echo "energy_rt plot failed (non-fatal)"

# the anisotropy observable: per-sector front radius r_q(theta)
$PYTHON experiments/energy_polar.py $WAVE -o $OUT/$NAME-polar.png \
  || echo "energy_polar plot failed (non-fatal)"

# warp so the peak deflection reads as ~5% of the unit domain in the render.
# --surface: grid lines finer than pixels moire as tubes, fill them in (see ne18-27)
WARP=$($PYTHON -c "
import h5py, numpy as np
with h5py.File('$WAVE') as f: m = np.abs(f['VTKHDF/PointData/values'][:]).max()
print(f'{0.05/max(m, 1e-30):.3g}')")
experiments/netvis.py $WAVE --view iso --surface --bg white --color-rescale frame --warp-by values:0,1,2 --warp-scale $WARP \
  --color-by values:2 --show 0 -o $VID \
  || echo "netvis animation failed (non-fatal, e.g. no pvpython)"
experiments/netvis.py $WAVE --view top --surface --bg white --warp-by values:0,1,2 --warp-scale $WARP \
  --color-by values:2 --show 0 --frames 0.15 -o $IMG \
  || echo "netvis still failed (non-fatal, e.g. no pvpython)"
