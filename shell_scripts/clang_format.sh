#!/bin/bash


COL='\e[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR='\e[0m'     # text format 'standard'


## Names of files that will be formatted. (Placeholders may be used.)
FORMAT_FILES=("*.cxx" "*.hxx")


## Actual command that is run when script is ecexuted.
for file in ${FORMAT_FILES[*]}; do find .. -wholename ${file} -print | xargs clang-format -i; done
