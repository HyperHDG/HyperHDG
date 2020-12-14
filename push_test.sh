#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'


echo -e "${COL}${BOLD}Clean up ...${NOR}"
(set -x; make clean)

echo -e "\n${COL}${BOLD}Format the cxx and hxx files ...${NOR}"
(set -x; make format)

echo -e "\n${COL}${BOLD}Do the GitHub run tests ...${NOR}"
(set -x; make test_github)

echo -e "\n${COL}${BOLD}Do the GitHub run tests ...${NOR}"
(set -x; make test_all_compilers)

echo -e "\n${COL}${BOLD}Clean up ...${NOR}"
(set -x; make clean)

echo -e "\n${COL}${BOLD}Make doxygen ...${NOR}"
(set -x; make doxygen > doxygen/doxy_log.txt)

echo -e "\n${COL}${BOLD}Re-install and test all components ...${NOR}"
(set -x; make build)