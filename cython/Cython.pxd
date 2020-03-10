from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Cythonize.cxx" :
  cdef void hyCython(vector[string])
