#!/bin/bash
# Spatial convergence plot of the compatible (homogeneous-Dirichlet) wave4 arm from
# ne9-12-conv-x-hom.sh: e_rel vs h, one line per degree, reference slopes h^{p+1} = 2/4/6.
# Companion figure: e_dual (the dual pair (n, m), same h^{p+1} rates) as
# ne9-12-conv-dual-hom.{png,tex,csv}. Excluded there: the fine-mesh tails (deg 3: nx > 32,
# deg 5: nx > 16) where the endpoint dual recovery's dt-independent roundoff floor
# (~1e-8 at nx 32, growing ~h^-5; cf. the reverted endpoint-algebraic solve, 8ec51489)
# overtakes the h^{p+1} error.
# Writes $OUT.png (view: sxiv) and the standalone pgfplots pair $OUT.tex/$OUT.csv.
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-12-conv-x-hom}"
DUAL=${OUT_DIR}/ne9-12-conv-dual-hom
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, deg: .poly_deg, e_rel}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_rel \
        -w 'nx >= 4' \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylabel 'max rel. $L^2$ error' --ylim '1e-12,10' \
        --ref "2;32,16;5e-3|4;32,16;8e-7|6;32,16;6e-11" \
        --save "$OUT.png" --tikz "$OUT" --legend 'lower right' --nshow

jq -c '.Stdout | {nx, deg: .poly_deg, e_dual}' $OUT.json \
    | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_dual \
        -w 'nx >= 4 and not (deg == 3 and nx > 32) and not (deg == 5 and nx > 16)' \
        --log xy --xbase 2 --xlabel 'discretization size $h$\strut' \
        --ylabel 'max rel. dual $L^2$ error' --ylim '1e-12,10' \
        --ref "2;32,16;1e-3|4;32,16;2e-7|6;16,8;1e-9" \
        --save "$DUAL.png" --tikz "$DUAL" --legend 'lower right' --nshow
