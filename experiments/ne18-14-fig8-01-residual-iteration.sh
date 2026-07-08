#!/bin/bash
# usage: ne18-14-fig8-plot.sh <log.json> [<out.png>]
set -eo pipefail

NAME=$(basename -s .sh $0)
LOG=$1
IMG=${2:-$NAME.png}

### one flat record per iteration, normalizing the old -err_monitor key names
jq -c '.Stdout | .H as $H | .ksp_monitor[] | {it, enorm, H: $H}' $LOG \
  | python $(dirname $0)/plot.py -x it -y enorm -g H --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_{\text{E}}$' \
      --nshow --save $IMG
echo wrote $IMG
