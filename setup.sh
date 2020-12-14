#!/bin/bash


BOLD=$(tput bold)
COL='\e[0;97;41m'
NOR='\e[0m'


echo -e "${COL}${BOLD}Initialize git-based submodules ...${NOR}"
(set -x; git submodule update --init --recursive)

echo -e "\n${COL}${BOLD}Make doxygen ...${NOR}"
(set -x; cd doxygen; rm -rf html latex doxy_log.txt; doxygen Doxyfile > doxy_log.txt)

# echo -e "$\n{COL}${BOLD}Do the tests as if we were GitHub ...${NOR}"
# echo -e "${COL}Remove previous build and __pycache__ directories:${NOR}"
# (set -x; rm -rf build output */output __pycache__ */__pycache__)
# echo -e "${COL}Make new build directory:${NOR}"
# (set -x; mkdir -p build)
# echo -e "${COL}Configure:${NOR}"
# (set -x; cd build; cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="-DNOFILEOUT" \
#            -DNOPYTHONTESTS=True ..)
# echo -e "${COL}Build tests (C++):${NOR}"
# (set -x; cd build; make)
# echo -e "${COL}Run the tests (C++):${NOR}"
# (set -x; cd build; make test)

echo -e "\n${COL}${BOLD}Do the full testing and installing of components ...${NOR}"
echo -e "${COL}Remove previous build and __pycache__ directories:${NOR}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${COL}Make new build directory:${NOR}"
(set -x; mkdir -p build)
echo -e "${COL}Configure:${NOR}"
(set -x; cd build; cmake ..)
echo -e "${COL}Build the tests (C++):${NOR}"
(set -x; cd build; make)
echo -e "${COL}Run the tests (C++ & Python):${NOR}"
(set -x; cd build; make test)