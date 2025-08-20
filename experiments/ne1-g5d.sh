#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

jq -c -s 'reduce .[] as $item ({}; .["\($item.backend)_\($item.p)"] += $item) | .[]' \
    $OUTPUT-ne1-g2b.json $OUTPUT-ps-*-driver.log \
    | tee $OUTPUT-ne1-g5d.json \
    | experiments/ne1-graph.py \
        -x overlap -y cg_avg_time \
        --ylabel t \
        --group-by backend \
        --scatter \
        --title "average cg iteration time over number of subdomains" \
        --save $OUTPUT-ne1-g5d.png
