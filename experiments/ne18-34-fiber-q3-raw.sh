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
# coarse LU (not the default Cholesky): at p=31 the coarse matrix Z^T A Z is numerically
# indefinite in double precision (entries span ~1e12+ between soft welds and the x1e6
# rotation rigidities; narrow 4-edge hats) — MUMPS Cholesky factors it without complaint
# but the factor emits Inf on the first solve (both ne18-20 and the first ne18-21 attempt
# died with DIVERGED_NANORINF at it 0). LU pivoting absorbs it; verified on pde12.
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod -net2as_coarse_pc_type lu"
KSP="-ksp_monitor_yaml -ksp_monitor_yaml_enorm -ksp_rtol 1e-9 -ksp_max_it 400"
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

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/ne18-14-fig8-01-residual-iteration.sh experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# ne18-21's classic H sweep WITHOUT the x1e6 rotation-rigidity rescale: raw props put
# l_c = sqrt(EI/kGA) far below every H (deep bending regime, the (H/l_c)^2 stable-
# decomposition penalty) -- exactly what the -net2as_cb_degs 1,1,3,3,3,1 enrichment
# (commit 8fc400fa: C^1 bicubic plate pairs + rigid rotations in the coarse span) is
# meant to absorb. Two arms per p: plain q1 (expected to degrade badly; -ksp_max_it
# caps it) vs q3-enriched. Reference: ne18-21 rescaled-q1 = 33/51/63/69 its.
# Known risk (8fc400fa): the enrichment was INERT on the real slab (the ~100um z-tilt
# breaks planar plate pairs; z-squash x0.01 restored it, 248->58 its) -- this sweep
# documents the honest raw-props baseline either way. p = 3..31 = fig8's H^-1 = 4..32
# (p < sqrt(np) accepted deliberately for comparability with ne18-21/fig8).
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 2e-2"
SUB="--subdivide 128"

$PYTHON experiments/make_geo2.py -i domains/fiber-2026-05-20/net2/sca -o $FIBER2RAW $DIR $SUB | tee $LOGG

$PYTHON experiments/gortz_constants.py $FIBER2RAW --mu --cells 4 8 16 32 \
  --csv $OUT/gortz-fiber2raw.csv | tee $OUT/gortz-fiber2raw.txt

# raw-q1 arms may diverge or hit the iteration cap; parallel then exits nonzero, which
# must not kill the sweep before the log conversion and plot (arms are independent)
parallel --results $LOG --progress --bar -j1 \
  "echo H: 1/{=2 \$_+=1 =}; $MPIRUN -n $NP $BUILD/network -test {3} -domain $DOMAIN-{1}.geo.h5 $KSP $NET {4} -net2as_p {2}; echo domain: {1}; echo cb: {4}" \
  ::: fiber2raw ::: 3 7 15 31 ::: constant ::: "" "-net2as_cb_degs 1,1,3,3,3,1" \
  || echo "parallel exitcode: $? (some arms failed/capped, continuing)"

yq -io json -I0 '.Stdout |= from_yaml' $LOG

bash experiments/ne18-14-fig8-01-residual-iteration.sh $LOG $IMG \
  || echo "fig8 plot failed (non-fatal: the extra cb arm dimension may need manual splitting)"
