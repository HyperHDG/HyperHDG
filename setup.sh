#!/bin/bash

BOLD=$(tput bold)
NORMAL=$(tput sgr0)
COL='\033[1;33m'
NC='\033[0m' # No Color


echo -e "${COL}${BOLD}Initialize git-based submodules ...${NC}"
(set -x; git submodule update --init --recursive)

echo -e "${COL}${BOLD}\nDo the tests as if we were GitHub ...${NC}"
echo -e "${COL}Remove previous build and __pycache__ directories:${NC}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${COL}Make new build directory:${NC}"
(set -x; mkdir -p build)
echo -e "${COL}Do the cmake:${NC}"
(set -x; cd build; cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="-DNOFILEOUT" \
           -DNOPYTHONTESTS=True ..)
echo -e "${COL}Build the executables:${NC}"
(set -x; cd build; make)
echo -e "${COL}Conduct the tests:${NC}"
(set -x; cd build; make test)

echo -e "${COL}${BOLD}\nDo the full testing and installing of components ...${NC}"
echo -e "${COL}Remove previous build and __pycache__ directories:${NC}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${COL}Make new build directory:${NC}"
(set -x; mkdir -p build)
echo -e "${COL}Do the cmake:${NC}"
(set -x; cd build; cmake ..)
echo -e "${COL}Build the executables:${NC}"
(set -x; cd build; make)
echo -e "${COL}Conduct the tests:${NC}"
(set -x; cd build; make test)