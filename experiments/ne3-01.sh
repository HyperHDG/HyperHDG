#!/bin/bash

mkdir -p output
parallel --progress --bar \
    'python experiments/heat.py -i {1} -n {2} -t {=3 $_=2**-$_ =} | tee output/ne3-01-i{1}-n{2}-t{3}.yaml' \
    ::: $(seq 6) ::: $(seq 6 13) ::: 0 1
yq -cr "." output/ne3-01*.yaml \
    | experiments/plot.py -f json -x iteration -y error -g timesteps \
        --where "theta == .5" --log xy --trans "2.**-x,y" --xlabel h --xbase 2 \
        --ref "4;1,2;1e-4"
