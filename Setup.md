# How to setup HyperHDG?

There are two basic types of setting up HyperHDG at the moment: a regular installation and a Docker
container based installation. The following two paragraphs describe how to regularly install
HyperHDG, while the last paragraph deals with the Docker based installation.


## Regular installation

### Install required packages

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


### Obtain and install HyperHDG

To obtain HyperHDG, enter the directory you want to clone HyperHDG into and do the following steps:

1. Clone the repository and give it the name `your_name`:

   1. To clone this repository with https use

          $ git clone https://github.com/AndreasRupp/HyperHDG.git your_name

   2. To clone this repository with ssh use

          $ git clone git@github.com:AndreasRupp/HyperHDG.git your_name

   3. Download and extract the `.zip` file from https://github.com/AndreasRupp/HyperHDG.git

2. Enter directory using

       $ cd your_name

3. Execute the script `setup.sh` to install HyperHDG by

       $ ./shell_scripts/setup.sh

   If the default compiler of your system does not support `C-++20`, please obtain a compiler that
   does so and run

       $ CXX=compiler_name ./shell_scripts/setup.sh

   instead. This configures the default compiler of HyperHDG to be `compiler_name`.


With all these steps done and all tests of `setup.sh` passed, HyperHDG is ready to be used.



## Docker based installation

For the Docker based installation, you will need to have [Docker](https://www.docker.com/) installed
on your computer. Afterwards, obtain HyperHDG (see the first point in the above paragraph) and enter
the directory of HyperHDG. If you have `git` installed, run `make submodules`. Otherwise, obtain the
[Dockerfile](https://github.com/HyperHDG/docker) and copy it to `submodules/docker`.

To build the Docker image from the Dockerfile run

    make docker_build

in HyperHDG's main directory. Afterwards, you can use the `run` command illustrated in the README of
the [docker page](https://github.com/HyperHDG/docker) with the `<tag>` set to `hyperhdg_docker`. For
the usage of the Docker notebook, please also refer to the README.