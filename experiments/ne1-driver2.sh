#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
DRIVER="experiments/driver.py"
DOM="domains/fiber_network_14871.geo.bin.zstd"
BACKENDS="kahip naiveH"
LOGLEVEL="warning"
PFLAGS="--progress --eta"

echo "pre-assembling lhs matrix"
$DRIVER $DOM "" --mat-only --mat=$OUTPUT-mat.npz --log-level=warning
echo "driver: varying p"
parallel $PFLAGS $DRIVER $DOM $OUTPUT-pss-{1}-d{2}-{3}.dom.zstd -n '{=1 $_=2**$_ =}' --maxiter=300 --log-level=${LOGLEVEL} --log-file=$OUTPUT-pss-{1}-d{2}-{3}-driver.log --no-cg-progress --mat=$OUTPUT-mat.npz --log-data=delta:{2},backend:{3} \
  ::: $(seq 5) ::: 2 ::: ${BACKENDS}
