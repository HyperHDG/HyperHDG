#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'


echo -e "${COL}${BOLD}Initialize git-based submodules ...${NOR}"
(set -x; make submodules)

echo -e "\n${COL}${BOLD}Clean up ...${NOR}"
(set -x; make clean)

echo -e "\n${COL}${BOLD}Make doxygen ...${NOR}"
(set -x; make doxygen > doxygen/doxy_log.txt)

echo -e "\n${COL}${BOLD}Install and test all components ...${NOR}"
(set -x; make build)