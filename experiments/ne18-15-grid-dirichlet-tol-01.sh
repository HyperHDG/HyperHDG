#!/bin/bash
# usage: ne18-14-fig8-plot.sh <log.json> [<out.png>]
set -eo pipefail

NAME=$(basename -s .sh $0)
LOG=$1
IMG=${2:-$NAME.png}

echo "number of edges the dirichlet border spans, diameter of subdomains, number of iterations"
jq -rs 'sort_by(.Stdout.m)[] | .Stdout | [.m,.net2as.H,.iterations] | @tsv' $LOG | column -t
