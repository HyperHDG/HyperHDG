HyperHDG can be set up natively on Linux, macOS, and Windows, or within a Docker container. For all
types of setup, please obtain HyperHDG first. We recommend to do so using the program `git`, which
you might need to install first (see the operating system specific sections below).

HyperHDG is a header-only C++ library. Using it from C++ only requires a C++20 compiler, CMake, and
LAPACK. Using it from Python additionally requires Python 3 with NumPy, SciPy, and Cython, since the
Python interface compiles the requested C++ classes just-in-time via Cython.

The combinations of operating systems and compilers that are regularly checked to work are those
of the continuous integration, cf. the files in [.github/workflows](../tree/main/.github/workflows):

| Operating system                  | Compiler                    | LAPACK                  | Python          |
|:----------------------------------|:----------------------------|:------------------------|:----------------|
| Ubuntu 24.04                      | `g++-14`, `clang++-19`      | `liblapack-dev`         | system packages |
| macOS (latest, Apple Silicon)     | Apple Clang                 | Accelerate framework    | pip             |
| Windows (MSYS2 UCRT64)            | MinGW-w64 `g++`             | OpenBLAS / LAPACK       | MSYS2 packages  |

Other compilers supporting C++20 (at least GCC 10 or Clang 10) are likely to work, too. Microsoft's
MSVC compiler is *not* supported, since HyperHDG uses GNU extensions such as `__attribute__`.


# Obtain HyperHDG

To obtain HyperHDG, enter the directory you want to clone HyperHDG into. Then, clone the repository
including its submodules and give it the name `your_name` using one of the following ways:

- Clone the repository using `git` and `https`:

       $ git clone --recursive https://github.com/HyperHDG/HyperHDG.git your_name

- Clone the repository using `git` and `ssh`:

       $ git clone --recursive git@github.com:HyperHDG/HyperHDG.git your_name

If you have cloned without `--recursive`, run `git submodule update --init --recursive` within the
directory `your_name` (the script `setup.sh` described below also does this for you).


# Install required packages

## Linux (Ubuntu)

On Ubuntu 24.04 LTS the required packages can be installed using

    $ sudo apt-get install git cmake g++ libblas-dev liblapack-dev \
      python3-dev python3-numpy python3-scipy python3-matplotlib python3-pandas python3-joblib \
      cython3 doxygen graphviz

This command installs

- `git`, `cmake`, and `g++` to obtain HyperHDG and to control and conduct the compilation. `cmake`
  controls the software compilation process in a compiler independent form. Instead of `g++`, you
  can also use `clang++` (version 10 or newer).

- `libblas-dev` and `liblapack-dev` allow to use the LAPACK library within HyperHDG. LAPACK contains
  efficient solvers for dense linear equation systems which are used to solve the element-local
  systems of equations defined by the local hybrid discontinuous Galerkin (HDG) solvers.

- `python3-dev`, `python3-numpy`, `python3-scipy`, `python3-matplotlib`, `python3-pandas`, and
  `python3-joblib` procure the necessary Python functionalities if Python scripts are used to run
  HyperHDG.

- `cython3` facilitates Python language extensions to run C or C++ code within Python scripts. Thus,
  this package is necessary if HyperHDG is run in Python scripts. Some (older) distributions might
  require `cython3` to be replaced by `cython`.

- `doxygen` and `graphviz` (optional) provide the documentation system for all C++ and Python code.
  These packages are used to create a Doxygen documentation on your computer. This is useful to look
  up interfaces and properties of functions utilized in HyperHDG. An on-line version of the
  documentation can also be found under this [link](https://hyperhdg.github.io/auto_pages/doxygen).

Other distributions provide equivalent packages (e.g. `dnf install gcc-c++ cmake lapack-devel
python3-devel python3-numpy python3-scipy python3-Cython` on Fedora). Instead of the system's Python
packages, you can also use a virtual environment as described for macOS below. In that case,
`python3-dev` (the Python headers) is still needed.

To visualize the output of simulations, we recommend to install [ParaView](https://www.paraview.org).
Jupyter notebooks in the directory `jupyter` additionally require JupyterLab (e.g. `pip install
jupyterlab` within a virtual environment).


## macOS

No package manager is required on macOS. Compiler, `git`, and LAPACK (as part of Apple's Accelerate
framework) are provided by Apple's Command Line Tools.

1. Install the Command Line Tools (Apple Clang, `git`, the macOS SDK) by

       $ xcode-select --install

   and confirm the dialog. If the dialog reports that the software is not available, install them
   from the terminal instead (the second command asks for your administrator password):

       $ touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress
       $ softwareupdate --list        # shows the label, e.g. "Command Line Tools for Xcode 27.0-27.0"
       $ sudo softwareupdate -i "Command Line Tools for Xcode 27.0-27.0"

   A full installation of Xcode works as well.

2. Create a Python virtual environment and install CMake and the Python packages into it:

       $ python3 -m venv ~/.venvs/hyperhdg
       $ source ~/.venvs/hyperhdg/bin/activate
       $ pip install cmake cython numpy scipy matplotlib pandas joblib

   Any Python 3 interpreter can be used for this (e.g. the one of the Command Line Tools, one from
   [python.org](https://www.python.org), Homebrew, `uv`, `pyenv`, or conda). Activate the virtual
   environment (second command) in every terminal in which you want to use HyperHDG.

3. (Optional) If you use [Homebrew](https://brew.sh), you can alternatively install `cmake`, and
   install `doxygen` and `graphviz` for the documentation via

       $ brew install cmake doxygen graphviz


## Windows

There are two ways to use HyperHDG on Windows:

- **Windows Subsystem for Linux (WSL 2)**: Install WSL 2 with Ubuntu (`wsl --install` in a
  PowerShell with administrator rights) and follow the Linux instructions within the Ubuntu
  terminal. This is the most convenient way if you are used to Linux.

- **Native Windows via MSYS2**: Install [MSYS2](https://www.msys2.org) and open the terminal
  *MSYS2 UCRT64*. Within this terminal, install all required packages using

      $ pacman -S --needed git mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake \
        mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-lapack \
        mingw-w64-ucrt-x86_64-python mingw-w64-ucrt-x86_64-cython \
        mingw-w64-ucrt-x86_64-python-numpy mingw-w64-ucrt-x86_64-python-scipy \
        mingw-w64-ucrt-x86_64-python-matplotlib mingw-w64-ucrt-x86_64-python-pandas \
        mingw-w64-ucrt-x86_64-python-joblib

  and select Ninja as CMake's build system by

      $ export CMAKE_GENERATOR=Ninja

  (add this line to the file `~/.bashrc` to make it permanent). All further commands have to be
  executed in the *MSYS2 UCRT64* terminal. Use the Python of MSYS2 (`/ucrt64/bin/python3`) to run
  HyperHDG scripts: The just-in-time compiled Python modules are built by MinGW-w64's `g++` and
  are only compatible with this Python, not with the one from python.org or the Microsoft Store.


# Setup HyperHDG

1. Enter the directory using

       $ cd your_name

2. Execute the script `setup.sh` to setup HyperHDG by

       $ ./shell_scripts/setup.sh

   If the default compiler of your system does not support `C++20`, please obtain a compiler that
   does so and run

       $ CXX=compiler_name ./shell_scripts/setup.sh

   instead. This configures the default compiler of HyperHDG to be `compiler_name`.

3. Follow the instructions given by the script and select your choice of setup. The script can also
   be run non-interactively: `./shell_scripts/setup.sh -C -d` builds HyperHDG on your computer and
   skips the Docker setup (`-c` skips the build, `-D` builds the Docker container).

The script initializes the submodules, removes previous builds, builds and runs all tests within
the directory `build`, and generates the Doxygen documentation (if `doxygen` is installed). With all
these steps done and all tests passed, HyperHDG is ready to be used.

Alternatively, you can conduct these steps manually using

    $ mkdir -p build && cd build
    $ cmake ..
    $ cmake --build .
    $ ctest --output-on-failure

CMake uses the Python interpreter that is active (e.g. of the activated virtual environment) or the
first `python3` found in your `PATH`. To use another interpreter, add
`-DPython3_EXECUTABLE=/path/to/python3` to the `cmake` command.


## Run a first example

After the setup, you can, for example, run the convergence test of the diffusion equation by

    $ python3 reproducibles_python/diffusion_convergence_elliptic.py

from the main directory of HyperHDG. The first run of each configuration takes some seconds, since
the required C++ classes are compiled just-in-time (subsequent runs reuse the compiled modules in
`build/shared_objects`). The errors are printed and stored in `output/diffusion_convergence_elliptic.txt`.


## Troubleshooting

- *Python version is not CMAKE's Python version!* or *Wrong Python version!*: The Python that runs
  the script differs from the one CMake found during the setup. Activate the correct environment, or
  remove the directory `build` and rerun the setup with the Python you want to use.

- *File containing definition of ... has not been found!*: The submodules are missing. Run
  `git submodule update --init --recursive`.

- Linker errors during the just-in-time compilation: The linker settings are taken from
  `build/cmake_cython.cfg`, which is generated by CMake. Check that LAPACK has been found by CMake
  (see the output of `cmake ..`) and adapt the file if necessary.


# Docker based setup

For the Docker based setup, you will need to have [Docker](https://www.docker.com/) installed on
your computer. Afterwards, obtain HyperHDG (see the [corresponding paragraph](#obtain-hyperhdg)) and
enter its directory. Run the shell script

    $ CXX=compiler_name ./shell_scripts/setup.sh -c -D

to setup HyperHDG within a Docker container. Here, `compiler_name` is the name of some C++ compiler
that supports C++20 (e.g. `g++-14`). The compiler need not be installed on your system; it will be
installed within the Docker container. However, you might need `root` privileges to run the command.

Afterwards, you can use the `run` command illustrated in the README of the HyperHDG [docker page](
https://github.com/HyperHDG/docker) with the `<tag>` set to `hyperhdg_docker`. For the usage of the
Docker notebook, please also refer to the README.
