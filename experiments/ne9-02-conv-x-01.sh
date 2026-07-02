#!/bin/bash
set -x

: "${OUT:=${OUT_DIR:=output}/ne9-02-conv-x}"
mkdir -p $OUT_DIR

jq -c '.Stdout | {nx, e_rel, deg: .poly_deg}' $OUT.json \
  | experiments/plot.py -g deg -x nx --trans '1/x,y' -y e_rel --log xy \
      --save "$OUT.png,$OUT.pgf" --tikz "$OUT" \
      --xbase 2 --xlabel 'discretization size h' --ylabel 'maximum $L^2$ error' \
      --legend '' --ref "2;64,32;4e-3|3;64,32;3e-5|4;32,16;3e-6"
