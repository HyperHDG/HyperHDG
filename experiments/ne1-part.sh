#!/usr/bin/env sh

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
PART="build-release/experiments/partitioner"
GRAPH="experiments/ne1-graph.py"
DOM="domains/fiber_network_14871.geo.bin.zstd"
P="100"
D="2"
PS="40"
DS="15"
BACKENDS="kahip naive"
LOGLEVEL="debug"
PFLAGS="--progress"

#### remove old data (be careful)
echo "$(ls $OUTPUT*)"
read -p "remove these files? " resp
if [[ $resp == [yY] ]]; then
  rm -rf $(ls $OUTPUT*)
else
  exit 1
fi


#### generate partitions

parallel $PFLAGS "$PART $DOM {1} $OUTPUT-{1}-d{2}-{3} --delta={2} --backend={3} --log-level=${LOGLEVEL} --log-file=$OUTPUT-ps-{1}-d{2}-{3}.log --square=true" \
  ::: $(seq 2 2 $PS) ::: $D ::: ${BACKENDS}
parallel $PFLAGS "$PART $DOM {1} $OUTPUT-{1}-d{2}-{3} --delta={2} --backend={3} --log-level=${LOGLEVEL} --log-file=$OUTPUT-ds-{1}-d{2}-{3}.log --square=true" \
  ::: $(seq 2 2 $PS) ::: $D ::: ${BACKENDS}

#### time/p

A=$OUTPUT-kahip-runid-partitions.json
B=$OUTPUT-kahip-runid-time.json
C=$OUTPUT-kahip-time-partitions.json

cat $OUTPUT-ps-*-kahip.log | jq -c 'select(.message == "args") | {runid, partitions}' \
    > $A
cat $OUTPUT-ps-*-kahip.log | jq -c 'select(.message == "partitioner_backend") | {runid, time}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x partitions -y time --save $OUTPUT-time-partitions.png --title "partitioner backend runtime (t) over number of subdomains (p)" < $C
