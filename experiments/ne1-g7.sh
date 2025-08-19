#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat $OUTPUT-ds-*-driver.log | \
    jq -c | tee $OUTPUT-ne1-g7.json | \
    experiments/ne1-graph.py \
      -x delta -y cg_avg_time \
      --xlabel '$\delta$' \
      --ylabel '$t_\text{avg}$' \
      --group-by backend \
      --title 'average cg iteration time over overlap parameter $\delta$' \
      --save $OUTPUT-ne1-g7.png
