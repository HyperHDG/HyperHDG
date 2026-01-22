#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/$(basename "${0%.sh}")}" # set default if unset
mkdir -p $OUT_DIR

if [ -z $NOGEN ]; then
   parallel --shuf --progress --bar --results $OUT.json \
    'build/rel/experiments/wave -plot 0 -dim 2 -nx {=1 $_=2**$_ =} -nt {=2 $_=2**$_ =}' \
    ::: $(seq 8) ::: 7 10 13
   echo "gen exit: $?"
   cp $OUT.json{,.$(date +%s)}
fi
[ -z $NOPLOT ] && yq -I0 -o=json ".Stdout | from_yaml" $OUT.json \
    | experiments/plot.py -f json -x h -y e_abs -g nt \
        --log xy --xbase 2 --ref "4;.25,.5;1e-5" --save $OUT.png,$OUT.pgf $PLOTARGS
