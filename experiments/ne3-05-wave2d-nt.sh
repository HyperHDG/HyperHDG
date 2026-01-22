#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

if [ -z $NOGEN ]; then
   parallel --shuf --progress --bar --results $OUT.json \
     'build/rel/experiments/wave -plot 0 -dim 2 -nx {1} -nt {=2 $_=2**$_ =} -theta {3} -deg {4}' \
     ::: 16 ::: 64 256 1024 ::: 1 .5 ::: 6
   echo "gen exit: $?"
   cp $OUT.json{,$(date +%s)}
fi
[ -z $NOPLOT ] && yq -I0 -o=json ".Stdout | from_yaml" $OUT.json \
    | experiments/plot.py -f json -x dt -y e_abs -g theta \
        --log xy --xbase 2 --save $OUT.png,$OUT.pgf \
        --ref "2;.0078125,.015625;1e-5" $PLOTARGS
