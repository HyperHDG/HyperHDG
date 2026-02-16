#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

if [ -z $NOGEN ]; then
  parallel -j 1 --progress --bar --results $OUT.json \
    'mpirun -n 8 build/rel/experiments/network -plot 0 -domain domains/paper-{1}.geo.h5 -mat_cache output/{1}-mat.bin -pc_type jacobi -ksp_monitor_yaml -ksp_atol 1e-10' \
      ::: small big
  echo "gen exit: $?"
  cp $OUT.json{,.$(date +%s)}
fi

[ -z $NOPLOT ] && yq -I0 -o=json '.Stdout |= from_yaml | .V[0] as $domain | .Stdout.ksp_monitor[] | {"it": .it, "rnorm": .rnorm, "domain": $d}' output/ne4-02.json \
    | experiments/plot.py -f json -x it -y rnorm -g domain --save $OUT.png $PLOTARGS --log xy --marker ''
