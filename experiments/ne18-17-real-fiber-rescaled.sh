#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
FIBER2T=$DOMAIN-fiber2t.geo.h5
FIBER2D=$DOMAIN-fiber2d.geo.h5
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
IMG=$OUT/$NAME.png
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
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
GORTZ=$OUT/gortz.txt

export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/ne18-14-fig8-01-residual-iteration.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-16 (fiber2 at original scale, extent 8000, unit stiffnesses) degrades ~linearly in
# the subdomain count: with l_c = sqrt(EI/GA) = 1 every H = ext/(p+1) >> l_c sits deep in
# the bending-dominated regime, where interpolating u and r componentwise creates shear
# energy ~ GA H^2/EI times the bending energy, so the q1 coarse space has no H-uniform
# stable decomposition (verified on scaled unit grids: flat at extent 1, blow-up >= extent
# ~10). Rescaling the network to the unit bbox selects H <= 1/3 <= l_c, where the
# Timoshenko trace system behaves like componentwise diffusion and the Gortz theory
# applies. Same geometry with diffusion runs as control (dirichlet mask 1 vs 63).
DIRT="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63"
DIRD="--dirichlet xmin=1 xmax=1 ymin=1 ymax=1"
SUB="--subdivide 128"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2T --rescale-bbox $DIRT --dirichlet-tol 2e-2 $SUB --no-props | tee $LOGG
$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2D --rescale-bbox $DIRD --dirichlet-tol 2e-2 $SUB --no-props | tee -a $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2T --mu --cells 3 4 8 10 12 14 | tee $GORTZ

parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=3 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {2} -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {3}; echo domain: {1}" \
  ::: fiber2t fiber2d :::+ constant diffusion ::: 2 3 7 11 13
echo "parallel exitcode: $?"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG
