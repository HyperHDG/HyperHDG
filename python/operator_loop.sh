#!/bin/bash
# Convergence and timing sweeps of python/operator_loop.py as independent processes.
#
#   bash python/operator_loop.sh          # from the repo root
#
# GNU parallel runs one solve per (dim, degree, n, loop, method) and stores every record in
# $OUT/<sweep>.json (one JSON object per line, the YAML of the solve parsed into .Stdout); jq
# turns the records into the tables printed at the end. Compiles run in parallel, the timing sweep
# runs serially so the numbers do not disturb each other.
set -eo pipefail

NAME=$(basename -s .sh "$0")
OUT=output/$NAME.$(date +%s)
PY="${PYTHON:-python} python/operator_loop.py"
mkdir -p "$OUT"
ln -sfn "$(basename "$OUT")" output/$NAME
git rev-parse HEAD > "$OUT/rev"

echo "== compile (parallel)"
parallel -j "$(nproc)" "$PY --dim {1} --degree {2} --n 1 > /dev/null" ::: 1 2 3 ::: 1 2 3

echo "== convergence (parallel): both loops, assembled"
parallel --results "$OUT/convergence.json" -j "$(nproc)" \
  "$PY --dim {1} --degree {2} --n {3} --loop {4}" ::: 1 2 ::: 1 2 3 ::: 1 2 4 8 16 ::: elliptic operator
yq -io json -I0 '.Stdout |= from_yaml' "$OUT/convergence.json"

echo "== timings (serial): assembled and matrix-free, best of 3"
printf '%s\n' "2 1 32" "2 3 16" "3 1 8" "3 2 4" "3 3 3" > "$OUT/cases"
parallel --results "$OUT/timings-operator.json" -j1 --colsep ' ' \
  "$PY --dim {1} --degree {2} --n {3} --loop operator --method {4} --prepare {5} --repeat 3" \
  :::: "$OUT/cases" ::: assembled matrixfree ::: true false
parallel --results "$OUT/timings-elliptic.json" -j1 --colsep ' ' \
  "$PY --dim {1} --degree {2} --n {3} --loop elliptic --method {4} --prepare {5} --repeat 3" \
  :::: "$OUT/cases" ::: assembled matrixfree ::: true
cat "$OUT/timings-operator.json" "$OUT/timings-elliptic.json" > "$OUT/timings.json"
yq -io json -I0 '.Stdout |= from_yaml' "$OUT/timings.json"

echo; echo "convergence: error of both loops per (dim, degree, n)"
{ echo "dim degree n dofs elliptic operator reldiff"
  jq -r '.Stdout | [.opt.dim, .opt.degree, .opt.n, .dofs, .opt.loop, .error] | @tsv' "$OUT/convergence.json" \
    | sort -n -k1,1 -k2,2 -k3,3 \
    | awk '{ k = $1" "$2" "$3" "$4; e[k, $5] = $6; if (!(k in seen)) { seen[k] = 1; order[++m] = k } }
           END { for (i = 1; i <= m; ++i) { k = order[i]; a = e[k, "elliptic"]; b = e[k, "operator"];
                 printf "%s %.6e %.6e %.1e\n", k, a, b, (a - b < 0 ? b - a : a - b) / a } }'
} | column -t | tee "$OUT/convergence.txt"

echo; echo "timings [ms]: setup (matrix or prepare), solve, one product, errors"
{ echo "dim degree n dofs loop method prepare iters t_setup t_solve t_apply t_errors"
  jq -r '.Stdout | [.opt.dim, .opt.degree, .opt.n, .dofs, .opt.loop, .opt.method, .opt.prepare, .iterations,
                    (.t_setup * 1000 | round), (.t_solve * 1000 | round), (.t_apply * 1e6 | round / 1000), (.t_errors * 1000 | round)] | @tsv' \
    "$OUT/timings.json" | sort -n -k1,1 -k2,2 -k5,5 -k6,6
} | column -t | tee "$OUT/timings.txt"
echo "records: $OUT/{convergence,timings}.json"
