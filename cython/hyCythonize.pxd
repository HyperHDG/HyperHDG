# C++: Cythonize.cxx  (This line allows for just in time compile optimization!)

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<Cythonize.cxx>" :
  cdef string hyCythonize(vector[string], unsigned int, unsigned int)
