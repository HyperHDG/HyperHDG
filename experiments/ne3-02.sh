#!/bin/bash
set -euox pipefail

mkdir -p output
parallel --progress --bar --results output/ne3-02.json \
    'build/rel/experiments/wave -nx {=1 $_=2**$_ =} -nt {=2 $_=2**$_ =}' \
    ::: $(seq 10) ::: $(seq 7 13)
yq -I0 -o=json ".Stdout | from_yaml" output/ne3-02.json \
    | experiments/plot.py -f json -x h -y e_abs -g nt \
        --log xy --xbase 2 --ref "4;.25,.5;1e-3"
