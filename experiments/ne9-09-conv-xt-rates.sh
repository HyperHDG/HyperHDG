#!/bin/bash
# Observed convergence orders from an ne9-09-conv-xt.sh results json: for each degree the
# e_rel ratio between successive x+t doublings, as order = log2(e_prev/e).
# Expect ~2 / 4 / 6 for deg 1 / 3 / 5 (stage count coupled to degree in timowave.cxx).
#   usage: experiments/ne9-09-conv-xt-rates.sh output/ne9-09-conv-xt.json
set -e
f="${1:?usage: ne9-09-conv-xt-rates.sh <results.json>}"
yq -p=json -r '[.Stdout.poly_deg, .Stdout.nx, .Stdout.e_rel] | @tsv' "$f" \
  | sort -n -k1,1 -k2,2 \
  | awk 'BEGIN { print "deg\tnx\te_rel\torder" }
         { o = ($1 == deg && prev > 0) ? sprintf("%.2f", log(prev/$3)/log(2)) : "-";
           printf "%d\t%d\t%.3e\t%s\n", $1, $2, $3, o; deg = $1; prev = $3 }'
