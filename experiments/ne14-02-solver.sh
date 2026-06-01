#!/bin/bash
set -xeo pipefail

: ${BUILD:=build/openblas/experiments}
: ${OUTDIR:=output}
: ${INPUT:=$HOME//networks/morgan-2026-05-20/net1/sca}
NAME=$(basename -s .sh $0)
NOW=$(date +%s)
OUT=$OUTDIR/$NAME.$NOW
DOMAIN=$OUT/domain.geo.h5
TRACE=$OUT/trace.h5
STATIC=$OUT/static.vtkhdf
WAVE=$OUT/wave.vtkhdf
VID=$OUT/$NAME.avi
IMG=$OUT/$NAME.png
PERF=$OUT/perf.flamegraph
LOG=$OUT/log.yaml
RES1=$OUT/res1.json
RES2=$OUT/res2.json
NET="-net2as_p 1 -net2as_cb_type pu -net2as_print_local -mem_max -net2as_pc_factor_mat_solver_type"
export OMP_NUM_THREADS=1

mkdir -p $OUT
ln -sfn $NAME.$NOW $OUTDIR/$NAME
cmake --build --preset openblas
cp $BUILD/{network,timowave} experiments/{make_geo2,netvis}.py $OUT
git rev-parse HEAD > $OUT/rev
if ! git diff-index --quiet HEAD; then echo '-dirty' >> $OUT/rev; fi

. venv/bin/activate
. spack/mkl/.spack-env/view/setvars.sh

parallel --progress --bar --results $RES1 \
  "python experiments/make_geo2.py -i $INPUT --clamp-xy {1} -o $OUT/domain-{2}.geo.h5 --dirichlet xmax=68 xmin=63" ::: $(seq .1 .2 1) :::+ $(seq 5)
echo "exit: $?"
parallel --progress --bar --results $RES2 --colsep ' ' \
  "echo blas: {2}; build/{2}/experiments/network -domain {1} -comp 2 -strain .15 $NET {3}" ::: $OUT/domain-1.geo.h5 :::: - <<EOF
openblas mumps
openblas cholmod
mkl mumps
mkl cholmod
mkl mkl_pardiso
EOF
echo "exit: $?"
yq -i '.Stdout |= from_yaml' $RES2

