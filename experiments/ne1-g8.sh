#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat output/test-ds-*-driver.log | \
    jq -c | tee $OUTPUT-ne1-g8.json | \
    experiments/ne1-graph.py \
      -x delta -y cg_total_time \
      --xlabel '$\delta$' \
      --ylabel '$t_\text{tot}$' \
      --group-by backend \
      --title 'total cg time over overlap parameter $\delta$' \
      --save $OUTPUT-ne1-g8.png
