#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

#### bal/p

OUTPUT=$1
A=$OUTPUT-ne1-g2-A.json
B=$OUTPUT-ne1-g2-B.json
C=$OUTPUT-ne1-g2.json
GRAPH="experiments/ne1-graph.py"

cat $OUTPUT-ps-*.log | jq -c 'select(.message == "args") | {runid, p: .partitions, backend}' \
    > $A
cat $OUTPUT-ps-*.log | jq -c 'select(.message == "bal after") | {runid, bal}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x p -y bal --group-by=backend --save $OUTPUT-ne1-g2.png --title "balance of partition over number of subdomains" < $C
