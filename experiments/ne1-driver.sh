#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
DRIVER="experiments/driver.py"
DOM="domains/fiber_network_14871.geo.bin.zstd"
BACKENDS="kahip naive naiveH"
LOGLEVEL="warning"
PFLAGS="--progress --eta --jobs 1"

echo "pre-assembling lhs matrix"
$DRIVER $DOM "" --mat-only --mat=$OUTPUT-mat.npz --log-level=warning
echo "driver: varying p"
parallel $PFLAGS "$DRIVER $DOM $OUTPUT-ps-{1}-d{2}-{3}.dom.zstd -n {1} --maxiter 1000 --log-level=${LOGLEVEL} --log-file=$OUTPUT-ps-{1}-d{2}-{3}-driver.log --no-cg-progress --mat=$OUTPUT-mat.npz --log-data=delta:{2},backend:{3}" \
  ::: $(seq 2 2 20) ::: 2 ::: ${BACKENDS}
echo "partitioner: varying delta"
parallel $PFLAGS "$DRIVER $DOM $OUTPUT-ds-{1}-d{2}-{3}.dom.zstd -n {1} --maxiter 1000 --log-level=${LOGLEVEL} --log-file=$OUTPUT-ds-{1}-d{2}-{3}-driver.log --no-cg-progress --mat=$OUTPUT-mat.npz --log-data=delta:{2},backend:{3}" \
  ::: 10 ::: $(seq 1 1 10) ::: ${BACKENDS}
