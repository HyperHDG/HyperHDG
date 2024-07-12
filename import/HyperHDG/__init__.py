## \package HyperHDG
#
#  \brief   Python interface for HyperHDG library.
#
#  This package is the main ingredient to use the HyperHDG library within Python scripts. Its
#  function \c import compiles the necessary components just-in-time and provides the as classes
#  that can be used within Python. To this end, the Cython package is heavily used.
#
#  \authors   Guido Kanschat, Heidelberg University, 2021.
#  \authors   Andreas Rupp, Heidelberg University, 2021.

import os, sys

# Functions thata are used to compile C++ code.
from .import_cxx.config import config
from .import_cxx.include import include

try:
  import fiber_network
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../submodules/fiber_network.git")
  import fiber_network
