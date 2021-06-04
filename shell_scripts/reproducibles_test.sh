#!/bin/bash


COL='\e[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR='\e[0m'     # text format 'standard'

TIME=15s
FILES=$(dirname $(readlink -f "$0"))/../reproducibles_python/*.py


cd $(dirname $(readlink -f "$0"))/..; rm -f output/reproducibles_test.txt
for file in $FILES; do
  echo "\n\nTry $file ..." |& tee -a output/reproducibles_test.txt
  timeout $TIME python3 $file True |& tee -a output/reproducibles_test.txt
done

echo -e "\n${COL}Check whether tests have passed ...${NOR}"
python3 ./shell_scripts/check_reproducibles_test.py
