# C++: Cythonize.cxx  (This line allows for just in time compile optimization!)

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<cythonize.cxx>" :
  cdef string hyCythonize(vector[string], vector[string], unsigned int, unsigned int, unsigned int)
