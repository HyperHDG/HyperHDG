#!/bin/bash
# Plot the gortz.pdf figure-8 reproduction: per-iteration K-norm error of PCG with the
# net2as subspace decomposition preconditioner, one curve per coarse scale H = 1/Hinv,
# from the newline-delimited JSON log written by ne18-14-fig8-{grid,fiber}.sh
# (-ksp_monitor_yaml -ksp_monitor_yaml_enorm; also reads the old -err_monitor format).
# Also prints Table-2-style (average, worst) convergence rates tau_(l) = e_l / e_{l-1}
# for l >= 2 (gortz.pdf eq. 6.2).
# usage: ne18-14-fig8-plot.sh <log.json> [<out.png>]
set -eo pipefail

LOG=$1
IMG=${2:-fig8.png}

{ printf 'H\titers\tcond\ttau_avg\ttau_max\n'
  jq -r '(.ksp_monitor // .ksp_err_monitor) as $mon
    | [$mon[] | (.enorm // .err_K)] as $e
    | ([range(1; $e|length) | select($e[.] > $e[0]*1e-12) | $e[.]/$e[.-1]][1:]) as $r
    | ["1/\(.Hinv // .net2as_p)", ($mon | last | .it), (.cond // "nan"),
       (if $r|length > 0 then ($r|add/length*1000|round/1000), ($r|max*1000|round/1000)
        else "nan", "nan" end)]
    | @tsv' $LOG
} | column -t

# one flat record per iteration, normalizing the old -err_monitor key names
jq -c '{Hinv: (.Hinv // .net2as_p)} + ((.ksp_monitor // .ksp_err_monitor)[] | {it, enorm: (.enorm // .err_K)})' $LOG \
  | python $(dirname $0)/plot.py -x it -y enorm -g Hinv --log y --marker "" \
      --xlabel "iteration" --ylabel '$\|u - u^{(\ell)}\|_K$' \
      --nshow --save $IMG
echo saved $IMG
