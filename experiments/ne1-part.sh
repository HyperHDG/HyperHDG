#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
PART="build-release/experiments/partitioner"
DOM="domains/fiber_network_14871.geo.bin.zstd"
BACKENDS="kahip naive naiveH"
LOGLEVEL="debug"
PFLAGS="--progress --eta --jobs 1"

echo "partitioner: varying p"
parallel $PFLAGS "$PART $DOM {1} $OUTPUT-ps-{1}-d{2}-{3} --delta={2} --backend={3} --kahip-mode=0 --log-level=${LOGLEVEL} --log-file=$OUTPUT-ps-{1}-d{2}-{3}.log --square=true" \
  ::: $(seq 2 2 20) ::: 2 ::: ${BACKENDS}

echo "partitioner: varying delta"
parallel $PFLAGS "$PART $DOM {1} $OUTPUT-ds-{1}-d{2}-{3} --delta={2} --backend={3} --kahip-mode=0 --log-level=${LOGLEVEL} --log-file=$OUTPUT-ds-{1}-d{2}-{3}.log --square=true" \
  ::: 10 ::: $(seq 0 1 10) ::: ${BACKENDS}
