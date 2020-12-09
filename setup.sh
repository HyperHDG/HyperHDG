#!/bin/bash


BOLD=$(tput bold)
COL='\033[0;33m'
NOR='\033[0m'


echo -e "${COL}${BOLD}Initialize git-based submodules ...${NOR}"
(set -x; git submodule update --init --recursive)

echo -e "${COL}${BOLD}\nDo the tests as if we were GitHub ...${NOR}"
echo -e "${COL}Remove previous build and __pycache__ directories:${NOR}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${COL}Make new build directory:${NOR}"
(set -x; mkdir -p build)
echo -e "${COL}Do the cmake:${NOR}"
(set -x; cd build; cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="-DNOFILEOUT" \
           -DNOPYTHONTESTS=True ..)
echo -e "${COL}Build the executables:${NOR}"
(set -x; cd build; make)
echo -e "${COL}Conduct the tests:${NOR}"
(set -x; cd build; make test)

echo -e "${COL}${BOLD}\nDo the full testing and installing of components ...${NOR}"
echo -e "${COL}Remove previous build and __pycache__ directories:${NOR}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${COL}Make new build directory:${NOR}"
(set -x; mkdir -p build)
echo -e "${COL}Do the cmake:${NOR}"
(set -x; cd build; cmake ..)
echo -e "${COL}Build the executables:${NOR}"
(set -x; cd build; make)
echo -e "${COL}Conduct the tests:${NOR}"
(set -x; cd build; make test)