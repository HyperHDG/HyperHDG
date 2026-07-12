#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOM=$OUT/net3raw.geo.h5
WAVE=$OUT/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png

# FINAL rung of the ripple ladder (ne18-27 grid, ne18-28 mikado, ne18-30 fiber2): the
# 100M demo -- drumhead tap on net3 (40mm x 40mm, 40.9M nodes / ~70M edges raw, same
# kg-um-s units and material as net2, so the SAME c_shear = 1.722e9 um/s).
# The tap scales with the domain (5x net2): std_x = 500 um (2.5% of the half-extent,
# same relative size as the other rungs), std_t = std_x/c ~ 2.5e-7, tap = 5 std_t,
# dt = 1.25e-7 <= std_t/2.  center->rim = 20 mm -> arrival at 1.25e-6 + 1.16e-5;
# T = 1.2e-5 (NT=96): front at ~18.5 mm, just before the clamped rim.
# OUTPUT SIZING at ~10^8 edges: per-cell energies would be ~5 GB/frame -> plot_energy
# defaults OFF here (energy budget / E(r,t) for net3 wait for the in-situ histograms);
# values disp at stride 2 -> 49 frames x ~1.7 GB.  Memory: trace system ~5e8 dofs,
# matrix estimate ~350 GB -- pde11 (503G) is borderline, pde12 (1T) preferred; run a
# smoke (NT=4 STRIDE=1) first and watch RSS before committing to the full run.
: ${NT:=96}
: ${DT:=1.25e-7}
T=$(python3 -c "print($NT*$DT)")
: ${THETA:=0.5}
: ${DEG:=3}
: ${STD_X:=5e2}
: ${STD_T:=2.5e-7}
: ${TAP:=1.25e-6}
: ${ENERGY:=5e-4}
: ${STRIDE:=2}
: ${PLOT_VALUES:=disp}
: ${PLOT_PROPS:=none}
: ${PLOT_ENERGY:=0}

NET="-pc_type net2as -net2as_p 8 -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
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
# NFS collective-write fix, see ne18-27 (verified 2026-07-11); harmless on local disks
export OMPI_MCA_fs_ufs_lock_algorithm=1

# MPI-parallel rendering when a standalone ParaView is available (pde11: ~/opt;
# benchmarked 2026-07-11: 32 ranks render the mikado video in 1:30 vs 6:56 serial)
PVDIR=${PVDIR:-$HOME/opt/paraview/ParaView-5.13.3-osmesa}
if [ -x "$PVDIR/bin/pvbatch" ]; then
  NETVIS="env LD_LIBRARY_PATH=$PVDIR/lib/mesa LP_NUM_THREADS=4 VTK_SMP_MAX_THREADS=4 \
    $PVDIR/bin/mpiexec -np ${RENDER_NP:-32} $PVDIR/bin/pvbatch experiments/netvis.py"
else
  NETVIS=experiments/netvis.py
fi

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/netvis.py experiments/energy.py experiments/energy_rt.py \
   experiments/energy_polar.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-21/22 canonical build, full domain (no --clamp-xy)
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"
REG="--rescale-props 1,1,1,1,1e6,1e6,1e6,1,1,1,1,1,1,1,1,1,1"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net3/sca -o $DOM \
  $DIR $SUB $REG | tee $LOGG

$PYTHON experiments/gortz_constants.py $DOM --mu --cells 4 8 16 \
  --csv $OUT/gortz-net3raw.csv | tee $OUT/gortz-net3raw.txt

stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test drumhead -domain $DOM -deg $DEG \
  -theta $THETA -nt $NT -T $T \
  -std_x $STD_X -std_t $STD_T -tap_time $TAP -energy $ENERGY \
  -plot_stride $STRIDE -plot_values $PLOT_VALUES -plot_props $PLOT_PROPS \
  -plot_energy $PLOT_ENERGY \
  $KSP $NET -plot $WAVE -print_timestep | tee $LOG

$PYTHON experiments/energy.py $WAVE -o $NRG || echo "energy plot failed (non-fatal)"
$PYTHON experiments/energy_rt.py $WAVE -o $OUT/$NAME-rt.png \
  || echo "energy_rt plot failed (non-fatal)"
$PYTHON experiments/energy_polar.py $WAVE -o $OUT/$NAME-polar.png \
  || echo "energy_polar plot failed (non-fatal)"

# warp so the peak deflection reads as ~5% of the xy extent (um units here, so the
# 0.05*extent generalization matters); tubes ~ fiber width, sub-pixel at full view
WARP=$($PYTHON -c "
import h5py, numpy as np
with h5py.File('$WAVE') as f:
    m = np.abs(f['VTKHDF/PointData/values'][:]).max()
    p = f['VTKHDF/Points'][:, :2]
    ext = (p.max(0) - p.min(0)).max()
print(f'{0.05*ext/max(m, 1e-30):.3g}')")
$NETVIS $WAVE --view iso --bg white --color-rescale frame --warp-by values:0,1,2 \
  --warp-scale $WARP --color-by values:2 -r 8 --show 0 -o $VID \
  || echo "netvis animation failed (non-fatal, e.g. no pvpython)"
$NETVIS $WAVE --view top --bg white --warp-by values:0,1,2 --warp-scale $WARP \
  --color-by values:2 -r 8 --show 0 --frames $(python3 -c "print($T/2)") -o $IMG \
  || echo "netvis still failed (non-fatal, e.g. no pvpython)"
