#!/usr/bin/env bash
# Convergence of python/diffusion.py: errors and observed rates over polynomial degrees and
# refinement levels, i.e. what reproducibles_python/diffusion_convergence_elliptic.py does,
# from the shell and on a hypergraph read from file. Run from the repo root:
#
#   bash python/convergence.sh
#   bash python/convergence.sh domains/simplex_1_3.geo "1 2 3 4" "1 2 4 8"
#
# The refinement level divides the mesh size of every hyperedge, so the L2 error is expected
# to decay like h^(degree+1) and the rate column to approach degree + 1.
#
# The manufactured solution of TestParametersSinEllipt is u = sin(pi/2 x) with a right hand
# side taken from its second derivative along x, so it solves the network problem only on
# hyperedges parallel to the x axis (the default domain). On a graph with oblique edges --
# domains/simplex_1_2.geo, a fiber network -- the rate column plateaus instead.

set -euo pipefail

DOMAIN=${1:-domains/cross.geo}
DEGREES=${2:-"1 2 3"}
REFINEMENTS=${3:-"1 2 4 8 16 32"}

printf 'python/diffusion.py on %s\n\n' "$DOMAIN"
printf '%7s %11s %10s %14s %7s\n' degree refinement unknowns error rate

for degree in $DEGREES; do
  previous_error=""
  previous_refinement=""
  for refinement in $REFINEMENTS; do
    report=$(python python/diffusion.py "$DOMAIN" --degree "$degree" --refine "$refinement")
    unknowns=$(echo "$report" | sed -n 's/^unknowns: //p')
    error=$(echo "$report" | sed -n 's/^error: //p')

    # observed rate log(error_previous / error) / log(refinement / refinement_previous)
    rate=$(awk -v e="$error" -v p="$previous_error" -v r="$refinement" -v q="$previous_refinement" \
             'BEGIN { if (p == "" || e <= 0 || p <= 0) print "-";
                      else printf "%.2f", log(p / e) / log(r / q) }')
    printf '%7s %11s %10s %14s %7s\n' "$degree" "$refinement" "$unknowns" "$error" "$rate"

    previous_error=$error
    previous_refinement=$refinement
  done
  echo
done
