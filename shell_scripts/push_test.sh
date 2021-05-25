#!/bin/bash


COL='\e[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR='\e[0m'     # text format 'standard'

TEST_COMPILER=("clang++-10" "clang++-11" "g++-10")


test_compiler()
{
  rm -rf build domains/*.pts.geo \
    doxygen/html doxygen/latex doxygen/doxy_log.txt \
    .ipynb_checkpoints */.ipynb_checkpoints */*/.ipynb_checkpoints */*.nbconvert.* jupyter/*.py \
    __pycache__ */__pycache__ */*/__pycache__ &&
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=$1 && cmake --build . --config Debug &&
  ctest -C Debug && cd ..
}



cd $(dirname $(readlink -f "$0"))/..

rm -rf build domains/*.pts.geo \
  doxygen/html doxygen/latex doxygen/doxy_log.txt \
  .ipynb_checkpoints */.ipynb_checkpoints */*/.ipynb_checkpoints */*.nbconvert.* jupyter/*.py \
  __pycache__ */__pycache__ */*/__pycache__ |& tee output/push_test.txt

./shell_scripts/clang_format.sh |& tee -a output/push_test.txt

for comp in ${TEST_COMPILER[*]}; do test_compiler ${comp} |& tee -a output/push_test.txt; done

rm -rf build domains/*.pts.geo \
  doxygen/html doxygen/latex doxygen/doxy_log.txt \
  .ipynb_checkpoints */.ipynb_checkpoints */*/.ipynb_checkpoints */*.nbconvert.* jupyter/*.py \
  __pycache__ */__pycache__ */*/__pycache__ |& tee output/push_test.txt

cd doxygen && doxygen Doxyfile |& tee -a ../output/push_test.txt && cd .. &&

./shell_scripts/setup.sh -Cd |& tee -a output/push_test.txt

echo -e "\n${COL}Check whether tests have passed ...${NOR}"
python3 ./shell_scripts/check_push_test.py
