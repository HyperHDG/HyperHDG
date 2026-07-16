#!/bin/bash
# Spatial convergence plot of the compatible (homogeneous-Dirichlet) wave4 arm from
# ne9-12-conv-x-hom.sh: e_rel vs h, one line per degree, reference slopes h^{p+1} = 2/4/6.
# Writes $OUT.png (view: sxiv) and the standalone pgfplots pair $OUT.tex/$OUT.csv.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-12-conv-x-hom}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_rel}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_rel \
        --log xy --xbase 2 --xlabel 'discretization size $h$' \
        --ylabel 'max rel. $L^2$ error' \
        --ref "2;16,32;6e-2|4;16,32;2.5e-5|6;16,32;7e-9" \
        --save "$OUT.png" --tikz "$OUT" --nshow
