#!/bin/bash


COL='\e[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR='\e[0m'     # text format 'standard'


## Names of files that will be formatted. (Placeholders may be used.)
FORMAT_FILES=("*.cxx" "*.hxx")
FORMAT_FOLDERS=("reproducibles_c++", "tests_c++", "include")

## Do nothing if not started from main directory.
if [ "$(dirname $(readlink -f "$0"))" = "$PWD/shell_scripts" ]; then
  ## Actual command that is run when script is ecexuted.
  for foler in ${FORMAT_FOLDERS[*]}; do
    for file in ${FORMAT_FILES[*]}; do 
      find ${folder} -wholename ${file} -print | xargs clang-format -i;
    done;
  done
else
  echo "You using this script from another directory than HyperHDG's main. This is not allowed."
fi
