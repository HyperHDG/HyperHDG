#!/bin/bash
set -x

# ne12-02 rerun after the trace-error norm fix (edge-length-weighted skeleton
# norm, commit b85937ea): on the uniform cross2 grid the new relative e_trace
# equals the old one times sqrt(h)/sqrt(2), so the spatial trace error is now
# ~ {95, 6.9, 0.385}/sqrt(2) * h^(deg+2) (constants measured in the OLD norm
# for test 4 on cross2).  nt is recalibrated so the Crank-Nicolson error
# ~ 104*dt^2 stays at half the NEW spatial trace error:
#   nt = sqrt(2*sqrt(2)*104/const) * nx^((deg+2)/2),  floor 500,
# i.e. the deg-3 nx-64 job grows from 270k to ~906k steps.  Multi-rank MPI per
# job is not available here (ASCII .geo domains cannot be distributed,
# read_domain.hxx), but at ~6 ms/step the three deg-3 nx-64 jobs run
# concurrently under parallel: wall ~ 1.6 h.  Emitted most expensive first.

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

DATA_DIR=.
DOMAIN="cross2.geo"
DOMAIN="-domain $DATA_DIR/domains/$DOMAIN"
: "${PRESET:=rel}"
BIN_DIR=build/$PRESET/experiments
export OMP_NUM_THREADS=1

cmake --build --preset $PRESET --target timowave

joblist() {
  local nx=( 64     32     16    8    4   2)
  local nt1=(901    500    500   500  500 500)
  local nt2=(26800  6700   1700  500  500 500)
  local nt3=(906000 161000 28500 5100 900 500)
  for i in 0 1 2 3 4 5; do
    echo "${nx[i]} 3 ${nt3[i]}"
    echo "${nx[i]} 2 ${nt2[i]}"
    echo "${nx[i]} 1 ${nt1[i]}"
  done
}

parallel --progress --bar --results $OUT.json --colsep ' ' \
  "$BIN_DIR/timowave $DOMAIN -deg {2} -nx {1} -theta .5 -test wave4 -nt {3} -ksp_type preonly -pc_type lu -loc_lu_full -tau_s {4}" \
  :::: <(joblist) ::: 1 0 -1
echo "gen exit: $?"
yq -i '.Stdout |= from_yaml' $OUT.json
cp $OUT.json{,.$(date +%s)}
