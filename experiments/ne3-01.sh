#!/bin/bash
set -eox pipefail

mkdir -p output
[ -z $NOGEN ] && parallel --shuf --progress --bar --results output/ne3-01.json \
    "build/rel/experiments/heat -plot 0 -i {1} -ts {2}" \
    ::: $(seq 6) ::: $(seq 6 13)
[ -z $NOPLOT ] && yq -I0 -o=json ".Stdout | from_yaml" output/ne3-01.json \
    | experiments/plot.py -f json -x it -y e_abs -g timesteps \
        --where "theta == .5" --log xy --trans "2.**-x,y" --xlabel h --xbase 2 \
        --ref "4;2,1;1e-5" --save output/ne3-01.png,output/ne3-01.pgf $PLOTARGS
