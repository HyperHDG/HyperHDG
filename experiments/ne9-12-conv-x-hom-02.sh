#!/bin/bash
# Hybrid-variable (trace lambda) spatial convergence from the ne9-12-conv-x-hom.sh data:
# e_trace vs h per degree, guide slopes at p+2 = 3/5/7 -- the trace superconverges with
# order AT LEAST p+2 (measured: exactly p+2 for p = 1, 2; p+2 plus ~0.5..1 for p = 3, 5,
# robust across profiles/BCs/refinement families, see ne9-12-conv-x-hom-03/-04 and
# ne9-13-conv-trace-stiff; mechanism unresolved). Excluded points: nx 2 (all mesh nodes on
# zeros of sin(2*pi*s), relative normalization breaks down) and the nt-256 temporal-floor
# tails (deg 3: nx > 32, deg 5: nx > 16). Right panel of the manuscript's shared-axis row:
# y limits pinned, tick labels suppressed. Writes ne9-12-conv-trace-hom.{png,tex,csv}.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-12-conv-x-hom}"
TRC=${OUT_DIR}/ne9-12-conv-trace-hom
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_trace}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_trace \
        -w 'nx > 2 and not (deg == 3 and nx > 32) and not (deg == 5 and nx > 16)' \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylim '1e-12,10' --ytickoff \
        --ref "3;64,32;1.5e-4|5;32,16;1e-9|7;16,8;1.5e-11" \
        --save "$TRC.png" --tikz "$TRC" --legend 'upper right' --nshow
