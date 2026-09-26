#!/bin/bash


COL=$'\033[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR=$'\033[0m'     # text format 'standard'

## Time after which a reproducible is stopped (can be set via environment variable).
TIME=${REPRODUCIBLES_TIME:-15s}
MAIN_DIR="$(cd "$(dirname "$0")" && pwd)/.."
FILES="${MAIN_DIR}"/reproducibles_python/*.py


## Scripts that are skipped, and why.
#  - fiber_network_*: need large domain files (fiber_network_615452*) that are not in the repository.
#  - advection_hypergraph: none of the linear solvers converges (also on main as of 2026-09).
#    Remove it from this list once this has been fixed.
SKIP="fiber_network_diffusion.py fiber_network_elastic.py advection_hypergraph.py"


cd "${MAIN_DIR}" && mkdir -p output && rm -f output/reproducibles_test.txt
for file in $FILES; do
  if [[ " ${SKIP} " == *" $(basename "$file") "* ]]; then
    echo "${COL}Skipping $(basename "$file") (cf. SKIP in reproducibles_test.sh).${NOR}"
    continue
  fi
  printf "\n\nTry %s ...\n" "$file" 2>&1 | tee -a output/reproducibles_test.txt
  timeout $TIME python3 "$file" True 2>&1 | tee -a output/reproducibles_test.txt
done

echo -e "\n${COL}Check whether tests have passed ...${NOR}"
python3 ./shell_scripts/check_reproducibles_test.py
