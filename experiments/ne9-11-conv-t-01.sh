#!/bin/bash
# Temporal convergence plot of the compatible (homogeneous-Dirichlet) wave4 arm from
# ne9-11-conv-t-stages.sh: e_rel vs dt at fixed spatial discretization (deg 5, nx 256),
# one line per stage count, reference slopes dt^{2s} = 2/4/6.
# Writes $OUT.png (view: sxiv) and the standalone pgfplots pair $OUT.tex/$OUT.csv.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-11-conv-t-stages-hom}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nt, stages: .gauss_stages, e_rel}' $OUT.json \
    | experiments/plot.py -g stages -x nt --trans '1/x,y' -y e_rel \
        --log xy --xbase 2 --xlabel 'time step $\Delta t$\strut' \
        --ylabel 'max rel. $L^2$ error' --ylim '1e-11,10' \
        --ref "2;32,16;4e-3|4;32,16;1.5e-6|6;32,16;4e-10" \
        --save "$OUT.png" --tikz "$OUT" --legend 'lower right' --nshow
