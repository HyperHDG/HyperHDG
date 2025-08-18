#!/usr/bin/env sh

if [[ -z $1 ]]; then
  echo "ERROR: usage: script.sh <OUTPUT>"
  exit 1
fi

OUTPUT=$1
VIEW="sxiv"

# gen data
experiments/ne1-part.sh $OUTPUT
experiments/ne1-driver.sh $OUTPUT

# gen plots
for i in $(seq 8); do
  experiments/ne1-g$i.sh $OUTPUT
done

VIEW $OUTPUT*.png
