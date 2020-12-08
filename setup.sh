#!/bin/bash

BOLD=$(tput bold)
NORMAL=$(tput sgr0)
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${RED}${BOLD}!! Initialize git-based submodules:${NC}"
(set -x; git submodule update --init --recursive)
echo -e "${RED}${BOLD}!! Remove previous build and __pycache__ directories:${NC}"
(set -x; rm -rf build output */output __pycache__ */__pycache__)
echo -e "${RED}${BOLD}!! Make new build directory:${NC}"
(set -x; mkdir -p build)
echo -e "${RED}${BOLD}!! Do the cmake:${NC}"
(set -x; cd build; cmake ..)
echo -e "${RED}${BOLD}!! Build the executables:${NC}"
(set -x; cd build; make)
echo -e "${RED}${BOLD}!! Conduct the tests:${NC}"
(set -x; cd build; make test)