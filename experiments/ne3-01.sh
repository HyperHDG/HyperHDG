#!/bin/bash

mkdir -p output
parallel --progress --bar \
    'python experiments/heat.py -i {1} -n {2} -t {=3 $_=2**-$_ =} | tee output/ne3-01-i{1}-n{2}-t{3}.yaml' \
    ::: $(seq 6) ::: $(seq 6 13) ::: 0 1
yq -cr "." output/ne3-01*.yaml \
    | experiments/plot.py -f json -x iteration -y error -g theta,timesteps --log xy --save output/ne3-01.png
