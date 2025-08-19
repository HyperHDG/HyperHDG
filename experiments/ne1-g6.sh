#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat $OUTPUT-ds-*-driver.log | \
    jq -c | tee $OUTPUT-ne1-g6.json | \
    experiments/ne1-graph.py \
      -x delta -y cg_iters \
      --xlabel '$\delta$' \
      --ylabel iters \
      --group-by backend \
      --title 'number of cg iterations over overlap parameter $\delta$' \
      --save $OUTPUT-ne1-g6.png
