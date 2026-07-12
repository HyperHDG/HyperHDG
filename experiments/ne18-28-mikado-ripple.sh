#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOM=$OUT/mikado.geo.h5
WAVE=$OUT/wave.vtkhdf
LOG=$OUT/log.yaml
LOGG=$OUT/make_geo.log
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
NRG=$OUT/$NAME-energy.png

# second rung of the output/render ladder (ne18-27 = grid 500^2): the SAME drumhead tap
# (identical std_x/std_t/tap_time/energy/T/NT/stride) on a random mikado network of
# comparable size, so the energy plots (budget, E(r,t)) compare directly against the
# ordered grid -- coherent anisotropic front there vs scattering + diffuse coda here.
# mass 1000 at fiber length r = 0.05: stick density mass*r = 50 ~ 8.9x the percolation
# threshold, ~20k fibers, ~300k nodes in the giant component.  Mean segment length
# ~ r/(1+k) ~ 1e-2 >! std_x = 6e-3: the tap still covers a few segments.  No subdivide
# needed: segments <= r = 0.05 < 1/p = 0.14 (Goertz R0 bound).
: ${MASS:=1000}
: ${R:=0.05}
: ${SEED:=0}
: ${NT:=600}
: ${T:=0.6}
: ${THETA:=0.5}
: ${DEG:=2}
: ${STD_X:=6e-3}
: ${STD_T:=2e-3}
: ${TAP:=1e-2}
: ${ENERGY:=5e-4}
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
# NFS collective-write fix, see ne18-27 (verified 2026-07-11); harmless on local disks
export OMPI_MCA_fs_ufs_lock_algorithm=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target timowave
cp $BUILD/timowave experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/netvis.py experiments/energy.py experiments/energy_rt.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

$PYTHON experiments/make_geo2.py --mikado $MASS $R --seed $SEED \
  --dirichlet xmin=63 xmax=63 ymin=63 ymax=63 -o $DOM | tee $LOGG

$PYTHON experiments/gortz_constants.py $DOM --mu --cells 4 8 16 \
  --csv $OUT/gortz-mikado.csv | tee $OUT/gortz-mikado.txt

stdbuf -oL $MPIRUN -n $NP $BUILD/timowave -test drumhead -domain $DOM -deg $DEG \
  -theta $THETA -nt $NT -T $T \
  -std_x $STD_X -std_t $STD_T -tap_time $TAP -energy $ENERGY \
  -plot_stride $STRIDE -plot_values $PLOT_VALUES -plot_props $PLOT_PROPS \
  $KSP $NET -plot $WAVE -print_timestep | tee $LOG

$PYTHON experiments/energy.py $WAVE -o $NRG || echo "energy plot failed (non-fatal)"
$PYTHON experiments/energy_rt.py $WAVE -o $OUT/$NAME-rt.png \
  || echo "energy_rt plot failed (non-fatal)"

# warp so the peak deflection reads as ~5% of the unit domain; thin tubes (mean mikado
# segment ~ 1e-2, radius well below that) keep the disordered texture readable
WARP=$($PYTHON -c "
import h5py, numpy as np
with h5py.File('$WAVE') as f: m = np.abs(f['VTKHDF/PointData/values'][:]).max()
print(f'{0.05/max(m, 1e-30):.3g}')")
RAD=1e-3
experiments/netvis.py $WAVE --view iso --bg white --color-rescale frame --warp-by values:0,1,2 --warp-scale $WARP \
  --color-by values:2 -r $RAD --show 0 -o $VID \
  || echo "netvis animation failed (non-fatal, e.g. no pvpython)"
experiments/netvis.py $WAVE --view top --bg white --warp-by values:0,1,2 --warp-scale $WARP \
  --color-by values:2 -r $RAD --show 0 --frames 0.15 -o $IMG \
  || echo "netvis still failed (non-fatal, e.g. no pvpython)"
