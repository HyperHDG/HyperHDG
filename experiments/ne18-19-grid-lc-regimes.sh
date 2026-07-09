#!/bin/bash
set -xeo pipefail

: ${PRESET:=openblas}
: ${BUILD:=build/$PRESET/experiments}
: ${OUTDIR:=output}
: ${NP:=$(nproc)}
: ${MPIRUN:=spack/$PRESET/.spack-env/view/bin/mpirun}
[ -x "$MPIRUN" ] || MPIRUN=mpirun
# spack env python when usable (working matplotlib on the pde cluster), else system python
if [ -z "${PYTHON:-}" ]; then
  PYTHON=spack/$PRESET/.spack-env/view/bin/python
  $PYTHON -c 'import numpy, h5py, pandas, matplotlib' >/dev/null 2>&1 || PYTHON=python
fi
export PYTHON
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain
LOG=$OUT/log.json
LOGG=$OUT/make_geo.log
NET="-pc_type net2as -net2as_cb_type q1 -net2as_cb_trim -net2as_pc_factor_mat_solver_type cholmod"
KSP="-ksp_monitor_yaml -ksp_rtol 1e-9"

export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset $PRESET --target network
cp $BUILD/network experiments/make_geo2.py experiments/gortz_constants.py \
   experiments/plot.py $0 $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo dirty >> $OUT/rev; fi

# Minimal demonstration of the two Timoshenko regimes of the q1 two-level ASM on the
# homogeneous unit grid (unit stiffnesses => bending length l_c = sqrt(EI/GA) = 1):
#
#   s1    grid 129 on the unit square:      H = 1/(p+1) <= 1/4  <= l_c  (shear regime)
#   s8000 same grid, coordinates x 8000:    H = 8000/(p+1) >= 250 >> l_c (bending regime)
#
# Coordinate scaling with fixed unit props is equivalent (up to the change of variables
# r -> s r, which leaves the preconditioned spectrum invariant because both the node-based
# subdomain spaces and the componentwise coarse space are invariant under per-component
# diagonal scaling) to dividing the rotation-block stiffnesses GI, EI by s^2 on the unit
# grid, i.e. to moving the single dimensionless number H/l_c.  Expectation: in the shear
# regime the coarse space gives flat kappa (vs the one-level H^-2 reference arm); in the
# bending regime coarse == nocoarse == one-level H^-2, i.e. the componentwise coarse
# space is inert (stable-decomposition penalty (H/l_c)^2, cf. gortz_constants l_c check).
#
#   g17 / g17s (both at scale 100, edge length 6.25 vs subdivided to 0.78 < l_c = 1):
# graph refinement does NOT change the regime -- l_c is a property of the operator
# EI/GA, not of the mesh; subdividing edges below l_c only makes each condensed element
# shear-dominated locally, the assembled chain keeps the same soft bending compliance
# over the subdomain scale H = 100/(p+1) >> l_c.  Expectation: g17 == g17s.
DIR="--dirichlet xmin=63 xmax=63 ymin=63 ymax=63 --dirichlet-tol 1e-3"

$PYTHON experiments/make_geo2.py --grid 129 --no-props $DIR -o $DOMAIN-s1.geo.h5    >  $LOGG
$PYTHON experiments/make_geo2.py --grid 129 --no-props $DIR -o $DOMAIN-s8000.geo.h5 --scale 8000 >> $LOGG
$PYTHON experiments/make_geo2.py --grid 17  --no-props $DIR -o $DOMAIN-g17.geo.h5   --scale 100  >> $LOGG
$PYTHON experiments/make_geo2.py --grid 17  --no-props $DIR -o $DOMAIN-g17s.geo.h5  --scale 100 --subdivide 128 >> $LOGG

for d in s1 s8000 g17 g17s; do
  $PYTHON experiments/gortz_constants.py $DOMAIN-$d.geo.h5 --cells 2 4 8 16 | tee $OUT/gortz-$d.txt
done

# arm A/B: both scales, coarse vs one-level reference.  NB --results must end in
# .json to get a json-lines file (anything else is taken as a directory prefix).
LOGAB=$OUT/log-ab.json
LOGC=$OUT/log-c.json
parallel --results $LOGAB --progress --bar -j1 \
  "echo Hinv: {=4 \$_+=1 =}; echo domain: {1}; echo mode: {2}; $MPIRUN -n $NP $BUILD/network -test constant -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {4} {3}" \
  ::: s1 s8000 ::: coarse nocoarse :::+ '' '-net2as_nocoarse' ::: 3 7 15 31
# arm C: graph refinement does not change the regime
parallel --results $LOGC --progress --bar -j1 \
  "echo Hinv: {=2 \$_+=1 =}; echo domain: {1}; echo mode: coarse; $MPIRUN -n $NP $BUILD/network -test constant -domain $DOMAIN-{1}.geo.h5 $KSP $NET -net2as_p {2}" \
  ::: g17 g17s ::: 3 7
echo "parallel exitcode: $?"

cat $LOGAB $LOGC > $LOG && rm $LOGAB $LOGC
yq -io json -I0 '.Stdout |= from_yaml' $LOG

jq -r '.Stdout | [.domain, .mode, .Hinv, .iterations, .cond] | @tsv' $LOG \
  | sort | column -t | tee $OUT/summary.txt

jq -c '.Stdout | select(.domain == "s1" or .domain == "s8000")
       | {Hinv, cond, series: (.domain + "-" + .mode)}' $LOG \
  | $PYTHON experiments/plot.py -x Hinv -y cond -g series --log xy --xbase 2 \
      --ref "2;8,32;40" --xlabel 'H^{-1}' --ylabel '$\kappa$' \
      --nshow --save $OUT/$NAME-regimes.png --tikz $OUT/$NAME-regimes

jq -c '.Stdout | select(.domain == "g17" or .domain == "g17s")
       | {Hinv, cond, series: .domain}' $LOG \
  | $PYTHON experiments/plot.py -x Hinv -y cond -g series --log xy --xbase 2 \
      --xlabel 'H^{-1}' --ylabel '$\kappa$' \
      --nshow --save $OUT/$NAME-subdiv.png --tikz $OUT/$NAME-subdiv
