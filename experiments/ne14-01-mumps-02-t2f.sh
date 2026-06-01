#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-01-mumps}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {t_t2f, m: .net2as.global.nz_mat}' $OUT/res2.json \
    | experiments/plot.py -x m -y t_t2f --log xy \
        --save "$OUT/ne14-01-mumps-02-t2f.png" \
        --xlabel 'number of edges' --ylabel 'trace to flux assembly time' \
        --ref "1;5e6,1.2e7;8e0"
#        --figsize '2.5,2.5' \
#        --xlabel 'timestep size $\Delta t$' --ylabel 'maximum $L^2$ error' \
#        --legend ''

