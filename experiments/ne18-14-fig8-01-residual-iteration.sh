#!/bin/bash
# usage: ne18-14-fig8-plot.sh <log.json> [<out.png>]
set -eo pipefail

NAME=$(basename -s .sh $0)
LOG=$1
DIR=$(dirname $LOG)
IMG=${2:-$DIR/$NAME.png}

### one flat record per iteration, normalizing the old -err_monitor key names
jq -c '.Stdout | {H,domain} + (.ksp_monitor[] | {it, enorm})' $LOG \
  | python $(dirname $0)/plot.py -x it -y enorm -g H --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG --group0 domain --tikz $DIR/$NAME
echo wrote ${IMG%.png}_domain*.png
