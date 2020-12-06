#!/bin/bash

BOLD=$(tput bold)
NORMAL=$(tput sgr0)
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${RED}${BOLD}!! Initialize git-based submodules:${NC}"
echo -e "${RED}   git submodule update --init --recursive${NC}"
git submodule update --init --recursive
echo -e "${RED}${BOLD}!! Remove previous build and __pycache__ directories:${NC}"
echo -e "${RED}   rm -rf build __pycache__${NC}"
rm -rf build __pycache__
echo -e "${RED}${BOLD}!! Make new build directory:${NC}"
echo -e "${RED}   mkdir -p build${NC}"
mkdir -p build
echo -e "${RED}${BOLD}!! Enter build directoy:${NC}"
echo -e "${RED}   cd build${NC}"
cd build
echo -e "${RED}${BOLD}!! Do the cmake:${NC}"
echo -e "${RED}   cmake ..${NC}"
cmake ..
echo -e "${RED}${BOLD}!! Build the executables:${NC}"
echo -e "${RED}   make${NC}"
make
echo -e "${RED}${BOLD}!! Conduct the tests:${NC}"
echo -e "${RED}   test${NC}"
make test