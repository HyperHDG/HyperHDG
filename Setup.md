There are two basic types of setting up HyperHDG at the moment: a regular installation and a Docker
container based installation. For both types, please obtain HyperHDG first. We recommend to do so
using the program `git`, which you might need to install first.


# Obtain HyperHDG

To obtain HyperHDG, enter the directory you want to clone HyperHDG into. Then, clone the repository
and give it the name `your_name` using one of the following ways:

- Clone the repository using `git` and `https`:

       $ git clone https://github.com/HyperHDG/HyperHDG.git your_name

- Clone the repository using `git` and `ssh`:

       $ git clone git@github.com:HyperHDG/HyperHDG.git your_name


# Regular installation

## Install required packages

Before you start using HyperHDG, you need to install some packages. Having Ubuntu 20.04 LTS as
operating system this can be done using

    $ sudo apt-get install doxygen graphviz cmake python3-dev python3-numpy python3-scipy cython3 \
      libblas-dev liblapack-dev


This command installs

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


## Install HyperHDG

1. Enter directory using

       $ cd your_name

2. Execute the script `setup.sh` to install HyperHDG by

       $ ./shell_scripts/setup.sh

   If the default compiler of your system does not support `C-++20`, please obtain a compiler that
   does so and run

       $ CXX=compiler_name ./shell_scripts/setup.sh

   instead. This configures the default compiler of HyperHDG to be `compiler_name`.

3. Follow the instructions given by the script and select your choice of setup.


With all these steps done and all tests of `setup.sh` passed, HyperHDG is ready to be used.


# Docker based installation

For the Docker based installation, you will need to have [Docker](https://www.docker.com/) installed
on your computer. Afterwards, obtain HyperHDG (see the [corresponding paragraph](#obtain-hyperhdg))
and enter its directory. Run the shell script

    $ CXX=compiler_name ./shell_scripts/setup.sh

and follow its instructions to setup HyperHDG within a Docker container.
Here, `compiler_name` is the name of some C++ compiler that support C++20 .The compiler need not be
installed on your system; it will be installed within the Docker container. However, you need `root`
privileges to run the command.

Afterwards, you can use the `run` command illustrated in the README of the HyperHDG [docker page](
https://github.com/HyperHDG/docker) with the `<tag>` set to `hyperhdg_docker`. For the usage of the
Docker notebook, please also refer to the README.
