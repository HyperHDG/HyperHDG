#!/usr/bin/env bash
# Run the local-operator prototype's convergence ladders and benchmarks as independent processes.
# Compiles (one nanobind module per dimension and degree) and the convergence ladders run in
# parallel with GNU parallel; the timed benchmarks run serially afterwards so they do not disturb
# each other. Run from the repo root:
#
#   bash python/operator_loop.sh            # serial build of the prototype
#   bash python/operator_loop.sh --openmp   # prototype compiled with OpenMP (OMP_NUM_THREADS)
#
# Results: output/operator_loop_convergence.txt, output/operator_loop_bench[_ompN].txt

set -euo pipefail
OPT=("$@")
PY="python python/operator_loop.py"

CONV_CASES=$(printf '%s\n' "1 1" "1 2" "1 3" "2 1" "2 2" "2 3")
BENCH_CASES=$(printf '%s\n' "2 1" "2 3" "3 1" "3 2" "3 3")

echo "== compiling modules in parallel"
printf '%s\n' "$CONV_CASES" "$BENCH_CASES" | sort -u | \
  parallel --colsep ' ' -j "$(nproc)" "$PY compile --dim {1} --degree {2} ${OPT[*]} > /dev/null"

echo "== convergence ladders in parallel"
printf '%s\n' "$CONV_CASES" | \
  parallel --colsep ' ' -j "$(nproc)" "$PY convergence --dim {1} --degree {2} ${OPT[*]} > /dev/null"
{ head -1 output/operator_loop_convergence_d1_p1.txt
  printf '%s\n' "$CONV_CASES" | while read -r d p; do
    tail -n +2 "output/operator_loop_convergence_d${d}_p${p}.txt"; done
} > output/operator_loop_convergence.txt
cat output/operator_loop_convergence.txt

echo "== benchmarks (serial)"
SUFFIX=""
for o in "${OPT[@]}"; do [ "$o" = "--openmp" ] && SUFFIX="_omp${OMP_NUM_THREADS:-$(nproc)}"; done
printf '%s\n' "$BENCH_CASES" | while read -r d p; do
  $PY bench --dim "$d" --degree "$p" "${OPT[@]}" > /dev/null
done
{ head -1 "output/operator_loop_bench_d2_p1${SUFFIX}.txt"
  printf '%s\n' "$BENCH_CASES" | while read -r d p; do
    tail -n +2 "output/operator_loop_bench_d${d}_p${p}${SUFFIX}.txt"; done
} > "output/operator_loop_bench${SUFFIX}.txt"
cat "output/operator_loop_bench${SUFFIX}.txt"
