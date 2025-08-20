#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

#### bal/p

OUTPUT=$1
A=$OUTPUT-ne1-g3-A.json
B=$OUTPUT-ne1-g3-B.json
C=$OUTPUT-ne1-g3.json
GRAPH="experiments/ne1-graph.py"

cat $OUTPUT-ds-*.log | jq -c 'select(.message == "args") | {runid, delta, backend}' \
    > $A
cat $OUTPUT-ds-*.log | jq -c 'select(.message == "bal after") | {runid, bal}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x delta --xlabel '$\delta$' -y bal --group-by=backend --save $OUTPUT-ne1-g3.png --title "balance of partition over overlap parameter $\delta$" < $C
