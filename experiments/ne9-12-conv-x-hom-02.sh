#!/bin/bash
# Hybrid-variable (trace lambda) spatial convergence from the ne9-12-conv-x-hom.sh data:
# e_trace vs h per degree. The nodal trace superconverges past h^{p+1}. nx 2 is excluded:
# every mesh node sits on a zero of sin(2*pi*s) there, so the relative normalization of
# e_trace breaks down (n_trace ~ round-off). Right panel of the manuscript's shared-axis
# row: y limits pinned, tick labels suppressed. Writes ne9-12-conv-trace-hom.{png,tex,csv}.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-12-conv-x-hom}"
TRC=${OUT_DIR}/ne9-12-conv-trace-hom
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_trace}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_trace -w 'nx > 2' \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylim '1e-12,10' --ytickoff \
        --ref "2.5;64,32;1.5e-4|5.5;32,16;1e-9|7.5;16,8;8e-12" \
        --save "$TRC.png" --tikz "$TRC" --legend 'lower right' --nshow
