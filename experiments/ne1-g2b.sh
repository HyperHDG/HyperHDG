#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

#### bal/p

OUTPUT=$1
A=$OUTPUT-ne1-g2b-A.json
B=$OUTPUT-ne1-g2b-B.json
C=$OUTPUT-ne1-g2b.json
GRAPH="experiments/ne1-graph.py"

cat $OUTPUT-ps-*.log | jq -c 'select(.message == "args") | {runid, p: .partitions, backend}' \
    > $A
cat $OUTPUT-ps-*.log | jq -c 'select(.message == "total overlap after") | {runid, overlap}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x p -y overlap --group-by=backend --save $OUTPUT-ne1-g2b.png --title "total overlap of dd over number of subdomains" < $C
