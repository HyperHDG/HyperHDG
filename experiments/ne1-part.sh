#!/usr/bin/env sh

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi
OUTPUT=$1
PART="build-release/experiments/partitioner"
GRAPH="experiments/ne1-graph.py"
DOM="domains/fiber_network_615452.geo.bin.zstd"
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

#### time/p

A=$OUTPUT-runid-partitions.json
B=$OUTPUT-runid-time.json
C=$OUTPUT-time-partitions.json

cat $OUTPUT-ps-*.log | jq -c 'select(.message == "args") | {runid, p: .partitions, backend}' \
    > $A
cat $OUTPUT-ps-*.log | jq -c 'select(.message == "partitioner_backend") | {runid, t: .time}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x p -y t --group-by=backend --save $OUTPUT-time-partitions.png --title "partitioner backend runtime (in s) over number of subdomains" < $C
