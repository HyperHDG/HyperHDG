#!/bin/bash


COL=$'\033[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR=$'\033[0m'     # text format 'standard'


## Version of clang-format that is used by the continuous integration (.github/workflows/
#  clang_format.yml). Other versions might format the code differently. Install it, e.g., via
#  "python3 -m pip install clang-format==18.1.8" (within a virtual environment).
CLANG_FORMAT_VERSION="18.1.8"
CLANG_FORMAT="${CLANG_FORMAT:-clang-format}"

## Names of files that will be formatted. (Placeholders may be used.)
FORMAT_FILES=("*.cxx" "*.hxx")
FORMAT_FOLDERS=("reproducibles_c++" "tests_c++" "include" "examples")

## Do nothing if not started from main directory.
if [ "$(cd "$(dirname "$0")" && pwd)" = "$PWD/shell_scripts" ]; then
  if ! command -v "${CLANG_FORMAT}" > /dev/null 2>&1; then
    echo "${COL}'${CLANG_FORMAT}' has not been found. Install it via:"
    echo "python3 -m pip install clang-format==${CLANG_FORMAT_VERSION}${NOR}"
    exit 1
  fi
  if ! "${CLANG_FORMAT}" --version 2> /dev/null | grep -q "version ${CLANG_FORMAT_VERSION}"; then
    echo "${COL}Warning: The continuous integration uses clang-format ${CLANG_FORMAT_VERSION}, but"
    echo "'${CLANG_FORMAT}' is: $("${CLANG_FORMAT}" --version 2> /dev/null || echo "not found")."
    echo "Install it via: python3 -m pip install clang-format==${CLANG_FORMAT_VERSION}${NOR}"
  fi
  ## Actual command that is run when script is ecexuted.
  for folder in ${FORMAT_FOLDERS[*]}; do
    [ -d "${folder}" ] || continue
    for file in ${FORMAT_FILES[*]}; do
      find "${folder}" -wholename "${file}" -print | xargs "${CLANG_FORMAT}" -i
    done
  done
else
  echo "You using this script from another directory than HyperHDG's main. This is not allowed."
fi
