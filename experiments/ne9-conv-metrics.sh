#!/bin/bash
# Extract the golden convergence metrics from a parallel-results json produced by
# ne9-01-conv-t.sh / ne9-02-conv-x.sh, as a stable sorted table for diffing a
# refactor against its pre-change baseline (behaviour-preserving cleanup check).
#   usage: experiments/ne9-conv-metrics.sh output/ne9-01-conv-t.json
# Columns: deg theta nt nx | e_rel e_trace iterations  (sorted, so order-independent).
set -e
f="${1:?usage: ne9-conv-metrics.sh <results.json>}"
yq -p=json -r \
  '[.Stdout.poly_deg, .Stdout.theta // "-", .Stdout.nt, .Stdout.nx, .Stdout.e_rel, .Stdout.e_trace, .Stdout.iterations] | @tsv' \
  "$f" | sort
