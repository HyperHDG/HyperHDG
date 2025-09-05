#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
PART="build-release/experiments/partitioner"
DOM="domains/fiber_network_14871.geo.bin.zstd"
BACKENDS="kahip naiveH"
LOGLEVEL="debug"
PFLAGS="--progress --eta"

echo "partitioner: varying p"
parallel $PFLAGS $PART $DOM '{=1 $_=2**$_ =}' $OUTPUT-pss-{1}-d{2}-{3} --delta={2} --backend={3} --kahip-mode=0 --log-level=${LOGLEVEL} --log-file=$OUTPUT-pss-{1}-d{2}-{3}.log \
  ::: $(seq 5) ::: 2 ::: ${BACKENDS}
