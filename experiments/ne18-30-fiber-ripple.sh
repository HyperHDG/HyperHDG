#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOM=$OUT/fiber2raw.geo.h5
WAVE=$OUT/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png

# third rung of the ripple ladder (ne18-27 grid, ne18-28 mikado): the SAME drumhead
# tap physics on the FULL fiber2 net (8mm x 8mm, 1.64M nodes / 2.36M edges raw,
# kg-um-s units), ne18-21/22 canonical props: RAW, rotation rigidities x1e6 (does NOT
# enter the wave speed -- gortz: c_shear = sqrt(kGA/(m/l)) = 1.722e9 um/s uniformly,
# p2 = p98), subdivide 128, deg 3 (deg 1 is ~11x too soft on fiber props, ne18-22).
# timing: center->rim = 4000 um -> arrival at tap_time + 4000/c = 2.57e-6 s.
#   std_x = 100 um  (~1 flexural wavelength, 1.2% of domain; ne17-05 precedent)
#   std_t = std_x/c ~ 5e-8 s, tap_time = 5 std_t, dt = 2.4e-8 <= std_t/2
#   T = 2.4e-6 (NT=100): front at ~3.7mm, just before the clamped rim.
# smoke first: NT=20 STRIDE=1 (T=4.8e-7, front at ~400 um) to check frames + energy.
: ${NT:=100}
: ${DT:=2.4e-8}
T=$(python3 -c "print($NT*$DT)")
: ${THETA:=0.5}
: ${DEG:=3}
: ${STD_X:=1e2}
: ${STD_T:=5e-8}
: ${TAP:=2.5e-7}
: ${ENERGY:=5e-4}
: ${STRIDE:=2}
: ${PLOT_VALUES:=disp}
: ${PLOT_PROPS:=none}

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

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $DOM \
  $DIR $SUB $REG | tee $LOGG

$PYTHON experiments/gortz_constants.py $DOM --mu --cells 4 8 16 \
  --csv $OUT/gortz-fiber2raw.csv | tee $OUT/gortz-fiber2raw.txt

stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test drumhead -domain $DOM -deg $DEG \
  -theta $THETA -nt $NT -T $T \
  -std_x $STD_X -std_t $STD_T -tap_time $TAP -energy $ENERGY \
  -plot_stride $STRIDE -plot_values $PLOT_VALUES -plot_props $PLOT_PROPS \
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
