#!/bin/bash
# Observed convergence orders from an ne9-09/ne9-10/ne9-11 results json: within each group
# the e_rel ratio between successive doublings of the sweep variable, as order =
# log2(e_prev/e). Groups default to poly_deg (expect ~2/4/6 for deg 1/3/5 with the matched
# stage coupling); pass a third argument to group by another Stdout field, e.g. gauss_stages
# for the ne9-11 fixed-degree stage ladder.
#   usage: experiments/ne9-conv-rates.sh <results.json> <nx|nt> [group-field=poly_deg]
set -e
usage="usage: ne9-conv-rates.sh <results.json> <nx|nt> [group-field=poly_deg]"
f="${1:?$usage}"
case "${2:?$usage}" in
  nx) k=2 ;; nt) k=3 ;;
  *) echo "$usage" >&2; exit 1 ;;
esac
g="${3:-poly_deg}"
yq -p=json -oy -r "[.Stdout.$g, .Stdout.nx, .Stdout.nt, .Stdout.e_rel] | @tsv" "$f" \
  | sort -n -k1,1 -k${k},${k} \
  | awk -v g="$g" 'BEGIN { print g "\tnx\tnt\te_rel\torder" }
         { o = ($1 == grp && prev > 0) ? sprintf("%.2f", log(prev/$4)/log(2)) : "-";
           printf "%d\t%d\t%d\t%.3e\t%s\n", $1, $2, $3, $4, o; grp = $1; prev = $4 }'
