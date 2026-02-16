#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

if [ -z $NOGEN ]; then
  parallel -j 1 --progress --bar --results $OUT.json \
    'build/rel/experiments/network -plot 0 -net2as_print_local -domain domains/paper-small.geo.h5 -mat_cache output/small-mat.bin -net2as_p {1} -net2as_pc_factor_mat_solver_type mumps' \
      ::: 1 2 4 8 16
  echo "gen exit: $?"
  cp $OUT.json{,.$(date +%s)}
fi
[ -z $NOPLOT ] && yq -I0 -o=json '.Stdout | from_yaml | .net2as | .sz as $sz | .local.[] | {"n": .size, "t": .time, "p": $sz}' output/ne4-01.json \
    | experiments/plot.py -f json --scatter -x n -y t -g p --save $OUT-a.png $PLOTARGS --log xy --ref "1;5e2,2e3;6e-4"
[ -z $NOPLOT ] && yq -I0 -o=json '.Stdout | from_yaml | .net2as | .sz as $sz | .local.[] | {"n": .size, "fill": .fill, "p": $sz}' output/ne4-01.json \
    | experiments/plot.py -f json --scatter -x n -y fill -g p --save $OUT-b.png $PLOTARGS --log xy
