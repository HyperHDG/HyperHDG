#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'

TIME=15s
FILES=$(dirname $(readlink -f "$0"))/../reproducibles_python/*.py


(set -x; cd $(dirname $(readlink -f "$0"))/..; rm -f output/reproducibles_test.txt)
for file in $FILES
do
  echo -e "\n\n${COL}${BOLD}Try $file ...${NOR}" |& tee -a output/reproducibles_test.txt
  (set -x; cd $(dirname $(readlink -f "$0"))/..; \
    timeout $TIME python3 $file True |& tee -a output/reproducibles_test.txt)
done

echo -e "\n${COL}${BOLD}Check whether tests have passed ...${NOR}"
(set -x; cd $(dirname $(readlink -f "$0")); python3 check_reproducibles_test.py)
