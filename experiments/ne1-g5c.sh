#!/usr/bin/env bash

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1

cat $OUTPUT-ne1-g1c.json $OUTPUT-ps-*-driver.log \
    | jq -c -s 'reduce .[] as $item ({}; .["\($item.p)_\($item.backend)"] += $item) | .[] | {t: .t + .cg_total_time, p, backend}' \
    | experiments/ne1-graph.py \
      -x p -y t\
      --ylabel '$t_\text{cg+part}$' \
      --group-by backend \
      --title "total cg iteration time plus partitioner time over number of subdomains" \
      --save $OUTPUT-ne1-g5c.png
