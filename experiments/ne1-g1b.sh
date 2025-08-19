#!/usr/bin/env bash

#### time/p

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1
A=$OUTPUT-ne1-g1b-A.json
B=$OUTPUT-ne1-g1b-B.json
C=$OUTPUT-ne1-g1b.json
GRAPH="experiments/ne1-graph.py"

cat $OUTPUT-ps-*.log | jq -c 'select(.message == "args") | {runid, p: .partitions, backend}' \
    > $A
cat $OUTPUT-ps-*.log | jq -c 'select(.message == "make_domains_overlap") | {runid, t: .time}' \
    > $B
jq -c -s 'reduce .[] as $item ({}; .[$item.runid] += $item) | .[]' \
    $A $B | tee $C

$GRAPH -x p -y t --log=y --group-by=backend --save $OUTPUT-ne1-g1b.png --title "delta-overlap runtime (in s) over number of subdomains" < $C
