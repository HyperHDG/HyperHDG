#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

#### bal/p

OUTPUT=$1
A=$OUTPUT-ne1-g2c-A.json
B=$OUTPUT-ne1-g2c-B.json
C=$OUTPUT-ne1-g2c.json
GRAPH="experiments/ne1-graph.py"

cat $OUTPUT-ds-*.log | jq -c 'select(.message == "args") | {runid, d: .delta, backend}' \
    > $A
cat $OUTPUT-ds-*.log | jq -c 'select(.message == "total overlap after") | {runid, overlap}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x d --xlabel '$\delta$' -y overlap --log=y --group-by=backend --save $OUTPUT-ne1-g2c.png --title "total overlap of dd over overlap parameter $\delta$" < $C
