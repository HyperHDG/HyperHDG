#!/bin/bash
set -eox pipefail

mkdir -p output
OUT=output/$(basename "${0%.sh}")

if [ -z $NOGEN ]; then
   mv $OUT.json $OUT.json.1 || true
   parallel --shuf --progress --bar --results $OUT.json \
     'build/rel/experiments/wave -plot 0 -nx {1} -nt {=2 $_=2**$_ =} -theta {3}' \
     ::: 1024 ::: $(seq 15) ::: 1 .5
fi
[ -z $NOPLOT ] && yq -I0 -o=json ".Stdout | from_yaml" $OUT.json \
    | experiments/plot.py -f json -x dt -y e_abs -g theta \
        --log xy --xbase 2 --save $OUT.png,$OUT.pgf \
        --ref "2;.0625,.25;1e-4" $PLOTARGS
