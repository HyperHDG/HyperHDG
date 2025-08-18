#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat output/test-ps-*-driver.log | \
    jq -c | tee $OUTPUT-ne1-g5.json | \
    experiments/ne1-graph.py \
      -x p -y cg_avg_time \
      --ylabel t \
      --group-by backend \
      --title "average cg iteration time over number of subdomains" \
      --save $OUTPUT-ne1-g5.png
