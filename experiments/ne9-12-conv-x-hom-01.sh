#!/bin/bash
# Spatial convergence plot of the compatible (homogeneous-Dirichlet) wave4 arm from
# ne9-12-conv-x-hom.sh: e_rel vs h, one line per degree, reference slopes h^{p+1} = 2/4/6.
# Writes $OUT.png (view: sxiv) and the standalone pgfplots pair $OUT.tex/$OUT.csv.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-12-conv-x-hom}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_rel}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_rel \
        -w 'nx >= 4' \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylabel 'max rel. $L^2$ error' --ylim '1e-12,10' \
        --ref "2;32,16;5e-3|4;32,16;8e-7|6;32,16;6e-11" \
        --save "$OUT.png" --tikz "$OUT" --legend 'lower right' --nshow
