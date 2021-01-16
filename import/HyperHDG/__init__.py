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

from .config import config
from .include import include
