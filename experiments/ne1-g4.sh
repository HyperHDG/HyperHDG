#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat output/test-ps-*-driver.log | \
    jq -c | tee $OUTPUT-ne1-g4.json | \
    experiments/ne1-graph.py \
      -x p -y cg_iters \
      --ylabel iters \
      --group-by backend \
      --title "number of cg iterations over number of subdomains" \
      --save $OUTPUT-ne1-g4.png
