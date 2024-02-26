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

# Functions thata are used to compile C++ code.
from .import_cxx.config import config
from .import_cxx.include import include

from .precond.gortz_hellman_malqvist_22 import gortz_hellman_malqvist_22