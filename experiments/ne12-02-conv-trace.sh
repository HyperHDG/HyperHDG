#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
BIN_DIR=build/rel/experiments
export OMP_NUM_THREADS=1

cmake --build --preset rel --target timowave

# nt is chosen per (deg, nx) such that the Crank-Nicolson error ~ 104*dt^2 stays at half the
# spatial trace error ~ {95, 6.9, 0.385} * h^(deg+1.5) (constants measured for test 4 on cross2),
# i.e. nt = sqrt(2*104/const) * nx^((deg+1.5)/2). A fixed temporal/spatial ratio keeps the EOC
# exact. Direct solve (preonly+lu) avoids the algebraic error floor of CG with rtol 1e-10, which
# polluted deg 3 at nx=64; it costs nothing since the matrix is factorized once.
# Emitted most expensive first so parallel schedules the long deg-3 jobs immediately.
joblist() {
  local nx=( 64     32    16    8    4   2)
  local nt1=(500    500   500   500  500 500)
  local nt2=(8000   2500  800   500  500 500)
  local nt3=(270000 57000 12000 2500 600 500)
  for i in 0 1 2 3 4 5; do
    echo "${nx[i]} 3 ${nt3[i]}"
    echo "${nx[i]} 2 ${nt2[i]}"
    echo "${nx[i]} 1 ${nt1[i]}"
  done
}

parallel --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg {2} -nx {1} -theta .5 -test wave4 -nt {3} -ksp_type preonly -pc_type lu -tau_s {4}" \
  :::: <(joblist) ::: 1 0 -1
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
