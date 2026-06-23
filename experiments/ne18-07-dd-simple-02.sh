#!/bin/bash
set -xeo pipefail

: ${OUTDIR:=output}
OUT=$OUTDIR/ne18-07-dd-simple
LOG=$OUT/log.json
IMG=$OUT/ne18-07-dd-simple-02-sdsize.png

# Distribution of net2as subdomain sizes vs p (smaller H = larger p).
# Baseline coarse only (nocoarse has identical subdomains; cb_trim drops the
# small boundary-ring subdomains and so changes the distribution).
# Split each run's net2as.local into nonzero sizes vs size-0 idle-rank slots:
# net2as keeps max(1, sz) slots per rank, so with fewer subdomains than ranks
# every idle rank emits one size-0 line -- those are placeholders, not empty
# network patches, so we count them (nidle) rather than plot them at 0.
jq -c 'select(.coarse=="")
       | (.net2as.local | map(.size)) as $s
       | {p: .net2as.p[0],
          sizes:  [$s[] | select(. >  0)],
          nidle: ([$s[] | select(. == 0)] | length)}' $LOG \
  | experiments/sdsize.py --save $IMG \
      --title "ne18-07 morgan net1: net2as subdomain sizes"
