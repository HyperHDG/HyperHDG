#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'


if [ $# -ne 1 ]; then
  echo "You have not set a compiler, please provide the compiler to be utilized."
  exit 1
fi

CXX=$1


echo -e "${COL}${BOLD}Clean up ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; make clean |& tee output/push_test.txt)

echo -e "\n${COL}${BOLD}Format the cxx and hxx files ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; make format |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Do the GitHub run tests ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; \
  make test_all_compilers |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Clean up ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; make clean |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Make doxygen ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; \
  make doxygen > doxygen/doxy_log.txt |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Re-install and test all components ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0"))/..; make build CXX=$CXX |& tee -a output/push_test.txt)

echo -e "\n${COL}${BOLD}Check whether tests have passed ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0")); python3 check_push_test.py)
