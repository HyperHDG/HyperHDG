#!/usr/bin/env bash
# Convergence of python/timowave.py on the wave4 manufactured solution (domains/cross2.geo):
# errors and observed rates of three sweeps. Run from the repo root:
#
#   bash python/timowave_convergence.sh            # all three sweeps
#   bash python/timowave_convergence.sh space tau  # a subset
#
#   space  spatial order, temporal error negligible: deg 1/3/5 with matched Gauss stages
#          s = (deg+1)/2, nt 256; expect e_rel ~ h^(deg+1) = 2/4/6. e_trace superconverges,
#          at least h^(deg+2). (Not nx 2: all its mesh nodes sit on zeros of the solution, which
#          makes the relative e_trace meaningless.)
#   time   temporal order, spatial error negligible: deg 5 at nx 32, s = 1/2/3 Gauss stages;
#          expect the classical Gauss order 2s = 2/4/6 (until the spatial floor ~1e-10).
#   tau    penalty scaling tau ~ h^s at deg 3 (s = 2, nt 256); expect e_rel ~ h^3 for
#          tau ~ h and h^4 for tau = 1 and (pre-asymptotically faster) for tau ~ 1/h.
#
# All sweeps use -wave4_px -pi/2: the standing wave then vanishes at every tip of the cross,
# the Dirichlet data is identically zero and Gauss collocation delivers its full order 2s. With
# time-dependent Dirichlet data (-wave4_px 0) the temporal order drops to the stiff order
# s+1 (odd s) / s (even s) -- a property of Gauss collocation, not of the discretization.
#
# The first run compiles each (degree, stages) instantiation once (a minute or two each).

set -euo pipefail

PX=-1.5707963267948966
SWEEPS=${*:-"space time tau"}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# sweep VARIABLE "VALUES" ARGS...: run timowave.py with ARGS plus -VARIABLE value for every value,
# print errors and observed rates w.r.t. that variable (it scales like 1/h resp. 1/dt)
sweep() {
  local var=$1 values=$2
  shift 2
  printf '\npython/timowave.py %s, varying -%s\n' "$*" "$var"
  printf '%6s %14s %6s %14s %6s %14s %6s\n' "$var" e_rel rate e_trace rate e_dual rate
  local prev="" prev_val=""
  for val in $values; do
    local report errs
    report=$(python python/timowave.py -wave4_px "$PX" "$@" "-$var" "$val")
    errs=$(echo "$report" | awk -F': ' '$1 == "e_rel" {r = $2} $1 == "e_trace" {t = $2}
                                        $1 == "e_dual" {d = $2} END {print r, t, d}')
    awk -v errs="$errs" -v prev="$prev" -v val="$val" -v pval="$prev_val" 'BEGIN {
      n = split(errs, e, " "); split(prev, p, " ")
      printf "%6s", val
      for (i = 1; i <= n; ++i) {
        if (pval == "" || e[i] <= 0 || p[i] <= 0) rate = "-"
        else rate = sprintf("%.2f", log(p[i] / e[i]) / log(val / pval))
        printf " %14s %6s", e[i], rate
      }
      printf "\n" }'
    prev=$errs
    prev_val=$val
  done
}

for s in $SWEEPS; do
  case $s in
    space)
      sweep nx "4 8 16 32" -deg 1 -nt 256
      sweep nx "4 8 16 32" -deg 3 -nt 256
      sweep nx "4 8 16" -deg 5 -nt 256
      ;;
    time)
      for stages in 1 2 3; do
        sweep nt "4 8 16 32" -deg 5 -nx 32 -stages "$stages"
      done
      ;;
    tau)
      for tau_s in 1 0 -1; do
        sweep nx "4 8 16 32" -deg 3 -nt 256 -tau_s "$tau_s"
      done
      ;;
    *)
      echo "unknown sweep '$s' (space, time, tau)" >&2
      exit 1
      ;;
  esac
done
