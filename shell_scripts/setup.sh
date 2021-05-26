#!/bin/bash


COL='\e[0;36m'  # text format 3 -> 'text color', 6 -> 'cyan'
NOR='\e[0m'     # text format 'standard'


## Function that sets up a clean version of the library.
#  The steps to achieve this goal are:
#  1. Find and enter the main directory of HyperHDG.
#  2. Initialize all submodules to be up to date with GitHub versions.
#  3. Delete all files that are not part of the library (including results, configurations, ...).
#  4. Build the library within the 'build' directory and tests its functions.
#  5. Build a version of the Doxygen of the version created by the aforementioned steps.
setup_library()
{
  cd $(dirname $(readlink -f "$0"))/.. &&

  git submodule update --init --recursive &&

  rm -rf build domains/*.pts.geo doxygen/html doxygen/latex doxygen/doxy_log.txt \
    .ipynb_checkpoints */.ipynb_checkpoints */*/.ipynb_checkpoints */*.nbconvert.* jupyter/*.py \
    output */output __pycache__ */__pycache__ */*/__pycache__ &&

  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Debug && cmake --build . --config Debug && ctest -C Debug && cd .. &&

  cd build/doxygen && doxygen Doxyfile && cd ../..
}


## Function that sets up a docker image named 'hyperhdg_docker' to include installation of HyperHDG.
setup_docker()
{
  docker build --build-arg INIT_COMMAND="apt-get install -y git doxygen graphviz cmake cython3 \
    libblas-dev liblapack-dev ipython3 $(CXX) && CXX=$(CXX) shell_scripts/setup.sh && cd build; \
    rm -r CMakeCache.txt CMakeFiles CTestTestfile.cmake Makefile Testing cmake_install.cmake \
    cython_files cython_log.txt examples shared_objects tests_c++ tests_python" \
    -f submodules/docker.git/Dockerfile -t hyperhdg_docker .
}


## Actual script that is run when file is executed.
#  It first checks whether HyperHDG should be installed on the computer. Afterwards, it checks
#  whether HyperHDG should be installed in a Docker container.
while getopts "cCdD" opt; do
  case $opt in
    (C) ynC=y;;
    (c) ynC=n;;
    (D) ynD=y;;
    (d) ynD=n;;
    (*) echo -e "${COL}Illegal option set.${NOR}" && exit 1;;
  esac
done

input_correct=false
while [ "$input_correct" = false ]; do
  if [ -z "$ynC" ]; then
    read -p "$(echo -e "${COL}Do you wish to install HyperHDG on your computer? [Yn] ${NOR}")" ynC
  fi
  if [ -z "$ynC" ]; then
    input_correct=true && setup_library
  else
    case $ynC in
      [Yy] | [Yy]es ) input_correct=true && setup_library;; 
      [Nn] | [Nn]o  ) input_correct=true && echo -e "${COL}Skipping installation ...${NOR}";;
      *             ) echo -e "${COL}Please answer yes or no.${NOR}";;
    esac
  fi
done

echo " "
echo -e "${COL}You might need to be root to install HyperHDG within a Docker container.${NOR}"
input_correct=false
while [ "$input_correct" = false ]; do
  if [ -z "$ynD" ]; then
    read -p "$(echo -e "${COL}Do you want to install HyperHDG in a Docker? [Yn] ${NOR}")" ynD
  fi
  if [ -z "$ynD" ]; then
    input_correct=true && setup_docker
  else
    case $ynD in
      [Yy] | [Yy]es ) input_correct=true && setup_docker;; 
      [Nn] | [Nn]o  ) input_correct=true && echo -e "${COL}Skipping installation ...${NOR}";;
      *             ) echo -e "${COL}Please answer yes or no.${NOR}";;
    esac
  fi
done
