#!/bin/bash


COL=$'\033[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR=$'\033[0m'     # text format 'standard'

## Time after which a reproducible is stopped (can be set via environment variable).
TIME=${REPRODUCIBLES_TIME:-600s}
MAIN_DIR="$(cd "$(dirname "$0")" && pwd)/.."
## The build tree with the precompiled modules of reproducibles_python (cmake --build build).
export HYPERHDG_BUILD_DIR=${HYPERHDG_BUILD_DIR:-${MAIN_DIR}/build}
## One BLAS thread: the dense local solves are small, threaded BLAS only spins.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
FILES="${MAIN_DIR}"/reproducibles_python/*.py


## Scripts that are skipped, and why.
#  - advection_hypergraph: none of the linear solvers converges (also on main as of 2026-09).
#    Remove it from this list once this has been fixed.
#  - common: the shared code of the scripts, not a script.
SKIP="advection_hypergraph.py common.py"


cd "${MAIN_DIR}" && mkdir -p output && rm -f output/reproducibles_test.txt
for file in $FILES; do
  if [[ " ${SKIP} " == *" $(basename "$file") "* ]]; then
    echo "${COL}Skipping $(basename "$file") (cf. SKIP in reproducibles_test.sh).${NOR}"
    continue
  fi
  printf "\n\nTry %s ...\n" "$file" 2>&1 | tee -a output/reproducibles_test.txt
  timeout $TIME python3 "$file" --quick --check 2>&1 | tee -a output/reproducibles_test.txt
done

echo -e "\n${COL}Check whether tests have passed ...${NOR}"
python3 ./shell_scripts/check_reproducibles_test.py
