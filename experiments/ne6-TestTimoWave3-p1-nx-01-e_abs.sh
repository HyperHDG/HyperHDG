: "${OUT:=${OUT_DIR:=output}/ne6-TestTimoWave3-p1-nx}" # set default if unset
mkdir -p $OUT_DIR

jq -cs 'map(.Stdout | {e_abs, nx}) | sort_by(.nx).[]' $OUT.json
jq -c '.Stdout | {e_abs, nx}' $OUT.json \
  | experiments/plot.py -x nx -y e_abs --log xy --ref '-1;1,3;3'
