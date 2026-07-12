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

# "ripple in the pond" (ne17-04) on an N^2 unit-square grid with unit synthetic props:
# first rung of the output/render ladder towards the 100M fiber demo (grid -> mikado ->
# fiber net2 -> net3).  The drumhead tap needs NO initial data (starts from rest, force
# only): f = amplitude * gauss(|x-center|; std_x) * gauss(t - tap_time; std_t) * e_z with
# amplitude = energy/(std_x^2 std_t), so narrowing the tap keeps the injected impulse
# fixed.  Scaling from ne17-04 (grid 100): h = 1/(N-1), std_x = 3h (point-like splash),
# std_t = std_x/3 (impulse over by tap_time + ~2 std_t), tap_time = 5 std_t, T = 0.6 ~
# one center->rim crossing (c ~ 1 on unit props) -- stop before reflections return.
# dt = T/NT = 1e-3 ~ h/2: the front crosses one cell in ~2 steps.
: ${N:=500}
: ${NT:=600}
: ${T:=0.6}
: ${THETA:=0.5}
: ${DEG:=2}
: ${STD_X:=6e-3}
: ${STD_T:=2e-3}
: ${TAP:=1e-2}
: ${ENERGY:=5e-4}
# new plot granularity flags (ne18-27ff): stride 10 -> 61 stored frames; values disp ->
# 3 of 18 components (renderers must then address the filtered layout: --warp-by
# values:0,1,2 --color-by values:2); per-cell energies kept at this scale for energy.py
# and for validating the upcoming in-situ histograms against a full-fidelity file.
: ${STRIDE:=10}
: ${PLOT_VALUES:=disp}
: ${PLOT_PROPS:=none}

NET="-pc_type net2as -net2as_p 8 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
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

# warp so the peak deflection reads as ~5% of the unit domain in the render.
# --surface: at N=500 the grid lines are denser than pixels and alias into moire as
# tubes; the filled Delaunay surface is the honest rendering at this density.
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
