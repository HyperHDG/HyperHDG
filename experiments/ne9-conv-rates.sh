#!/bin/bash
# Observed convergence orders from an ne9-09-conv-x.sh / ne9-10-conv-t.sh results json:
# for each degree the e_rel ratio between successive doublings of the sweep variable, as
# order = log2(e_prev/e). Expect ~2/4/6 for deg 1/3/5 (stage count coupled to the degree
# in timowave.cxx).
#   usage: experiments/ne9-conv-rates.sh <results.json> <nx|nt>
set -e
usage="usage: ne9-conv-rates.sh <results.json> <nx|nt>"
f="${1:?$usage}"
case "${2:?$usage}" in
  nx) k=2 ;; nt) k=3 ;;
  *) echo "$usage" >&2; exit 1 ;;
esac
yq -p=json -oy -r '[.Stdout.poly_deg, .Stdout.nx, .Stdout.nt, .Stdout.e_rel] | @tsv' "$f" \
  | sort -n -k1,1 -k${k},${k} \
  | awk 'BEGIN { print "deg\tnx\tnt\te_rel\torder" }
         { o = ($1 == deg && prev > 0) ? sprintf("%.2f", log(prev/$4)/log(2)) : "-";
           printf "%d\t%d\t%d\t%.3e\t%s\n", $1, $2, $3, $4, o; deg = $1; prev = $4 }'
