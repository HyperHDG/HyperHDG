#!/usr/bin/env sh

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
PART=${PART:-"build/experiments/partitioner"}
DOM=${DOM:-"domains/fiber_network_14871.geo.bin.zstd"}
P=${P:-"100"}
D=${D:-"2"}
PS=${PS:-"40"}
DS=${DS:-"15"}
BACKENDS=${BACKENDS:-"kahip naive"}
LOGLEVEL=${LOGLEVEL:-"debug"}

echo "output=$OUTPUT"
echo "part=$PART"
echo "dom=$DOM"
echo "pflags='$PFLAGS'"
echo "ps=(seq 2 2 $PS)"
echo "ds=(seq 2 2 $DS)"
echo "backends=${BACKENDS}"
echo "loglevel=${LOGLEVEL}"

parallel $PFLAGS "$PART $DOM {1} $OUTPUT-{1}-d{2}-{3} --delta={2} --backend={3} --log-level=${LOGLEVEL} --log-file=$OUTPUT-{1}-d{2}-{3}.log --square=true" ::: $(seq 2 2 $PS) ::: $D ::: ${BACKENDS}
#parallel $PFLAGS "$PART $DOM {1} $OUTPUT-{1}-d{2}-{3} --delta={2} --backend={3} --log-level=${LOGLEVEL} --log-file=$OUTPUT-{1}-d{2}-{3}.log --square=true" ::: $P ::: $(seq 2 2 $DS) ::: ${BACKENDS}
