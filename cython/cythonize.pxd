from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<cython/cythonize.cxx>" :
  cdef string cythonize(vector[string], vector[string], unsigned int, unsigned int, unsigned int)