#!/bin/bash
set -euox pipefail

mkdir -p output
parallel --progress --bar --results output/ne3-01.json \
    'build/rel/experiments/heat -i {1} -ts {2} -t {=3 $_=2**-$_ =}' \
    ::: $(seq 6) ::: $(seq 6 13) ::: 1
yq -I0 -o=json ".Stdout | from_yaml" output/ne3-01.json \
    | experiments/plot.py -f json -x it -y e_abs -g timesteps \
        --where "theta == .5" --log xy --trans "2.**-x,y" --xlabel h --xbase 2 \
        --ref "4;1,2;1e-4"
