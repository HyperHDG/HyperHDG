#!/bin/bash
# Plot the gortz.pdf figure-8 reproduction: per-iteration K-norm error of PCG with the
# net2as subspace decomposition preconditioner, one curve per coarse scale H = 1/Hinv,
# from the newline-delimited JSON log written by ne18-14-fig8-{grid,fiber}.sh
# (-ksp_monitor_yaml -ksp_monitor_yaml_enorm; also reads the old -err_monitor format).
# Also prints Table-2-style (average, worst) convergence rates tau_(l) = e_l / e_{l-1}
# for l >= 2 (gortz.pdf eq. 6.2).
# usage: ne18-14-fig8-plot.sh <log.json> [<out.png>]
set -eo pipefail

NAME=$(basename -s .sh $0)
LOG=$1
IMG=${2:-$NAME.png}

### one flat record per iteration, normalizing the old -err_monitor key names
jq -c '.Stdout | .H as $H | .ksp_monitor[] | {it, enorm, H: $H}' $LOG \
  | python $(dirname $0)/plot.py -x it -y enorm -g H --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_K$' \
      --nshow --save $IMG
echo wrote $IMG
