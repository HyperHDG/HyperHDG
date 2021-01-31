# How to setup HyperHDG?


## Install required packages

Before you start using HyperHDG, you need to install some packages. Having Ubuntu 20.04 LTS as
operating system this can be done using

    $ sudo apt-get install git doxygen graphviz cmake python3-dev python3-numpy python3-scipy \
    cython3 libblas-dev liblapack-dev


This command installs

- `git` is the revision control system, GitHub is based on. It is the recommended tool to obtain and
administrate HyperHDG.

- `doxygen` and `graphviz` provide the documentation system for all C++ and Python code. These
packages are used to create a Doxygen documentation on your computer. This is useful to look up 
interfaces and properties of functions utilized in HyperHDG. An on-line version of the documentation
can also be found under this [link](https://hyperhdg.github.io/auto_pages/doxygen).

- `cmake` is used to control the software compilation process in an compiler independent form.

- `python3-dev`, `python3-numpy`, and `python3-scipy` procure the necessary Python functionalities
if Python scripts are used to run HyperHDG.

- `cython3` facilitates Python language extensions to run C or C++ code within Python scripts. Thus,
this package is necessary if HyperHDG is run in Python scrips. However, some Ubuntu distributions
might require `cython3` to be replaced by `cython`.

- `libblas-dev` and `liblapack-dev` allow to use the LAPACK library within HyperHDG. LAPACK contains
efficient solvers for dense linear equation systems and which are used to solve the element-local
systems of equations defined by the local hybrid discontinuous Galerkin (HDG) solvers.


Compilation of the C++ code can be done using a compiler that can deal with the standard `C++20`. A
list of compilers that are regularly check to work can be found in the `Makefile`. Moreover, to
visualize the output of simulations, we recommend to install `ParaView`.


## Obtain and install HyperHDG

To obtain HyperHDG, enter the directory you want to clone HyperHDG into and do the following steps:

1. Clone the repository and give it the name `your_name`:

   1. To clone this repository with https use

          $ git clone https://github.com/AndreasRupp/HyperHDG.git your_name

   2. To clone this repository with ssh use

          $ git clone git@github.com:AndreasRupp/HyperHDG.git your_name

2. Enter directory using

       $ cd your_name

3. Execute the script `setup.sh` to install HyperHDG by

       $ ./shell_scripts/setup.sh

   If the default compiler of your system does not support `C-++20`, please obtain a compiler that
   does so and run

       $ CXX=compiler_name ./shell_scripts/setup.sh

   instead. This configures the default compiler of HyperHDG to be `compiler_name`.


With all these steps done and all tests of `setup.sh` passed, HyperHDG is ready to be used.
