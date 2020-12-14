#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'


echo -e "${COL}${BOLD}Clean up ...${NOR}"
(set -x; cd ..; make clean |& tee output/push_test.txt)

echo -e "\n${COL}${BOLD}Format the cxx and hxx files ...${NOR}"
(set -x; cd ..; make format |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Do the GitHub run tests ...${NOR}"
(set -x; cd ..; make test_github |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Do the GitHub run tests ...${NOR}"
(set -x; cd ..; make test_all_compilers |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Clean up ...${NOR}"
(set -x; cd ..; make clean |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Make doxygen ...${NOR}"
(set -x; cd ..; make doxygen > doxygen/doxy_log.txt |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Re-install and test all components ...${NOR}"
(set -x; cd ..; make build |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Check whether tests have passed ...${NOR}"
(set -x; cd ..; python3 chec_push_test.py)