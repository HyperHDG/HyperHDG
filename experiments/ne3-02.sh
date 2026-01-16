#!/bin/bash
set -eox pipefail

mkdir -p output
OUT=output/$(basename "${0%.sh}")

if [ -z $NOGEN ]; then
   mv $OUT.json $OUT.json.1 || true
   parallel --shuf --progress --bar --results $OUT.json \
    'build/rel/experiments/wave -plot 0 -nx {=1 $_=2**$_ =} -nt {=2 $_=2**$_ =}' \
    ::: $(seq 8) ::: $(seq 7 3 16)
fi
[ -z $NOPLOT ] && yq -I0 -o=json ".Stdout | from_yaml" $OUT.json \
    | experiments/plot.py -f json -x h -y e_abs -g nt \
        --log xy --xbase 2 --ref "4;.25,.5;1e-5" --save $OUT.png,$OUT.pgf $PLOTARGS
