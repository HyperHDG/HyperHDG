#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne14-01-mumps}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {mem_max, nz_fac: .net2as.local[0].nz_fac}' $OUT/res2.json \
    | experiments/plot.py -x mem_max -y nz_fac --log xy \
        --save "$OUT/ne14-01-mumps-03-mem.png" \
        --xlabel 'number of nz in factored matrix' --ylabel 'maximum memory usage' \
        --ref "1.2;1e9,3e9;5e7" \
#        --figsize '2.5,2.5' \
#        --legend ''

