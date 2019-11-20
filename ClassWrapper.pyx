# distutils: language = c++
# distutils: sources = DiffusionProblem.C

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector

# c++ interface to cython
cdef extern from "DiffusionProblem.h":
  cdef cppclass DiffusionProblem "DiffusionProblemRegular<1,3,1>":
        DiffusionProblem(vector[int]) except +
        void read_dirichlet_indices(vector[int]);
        vector[double] return_zero_vector()
        vector[double] matrix_vector_multiply(vector[double])
        int size_of_system()
        void plot_solution(vector[double])

# creating a cython wrapper class
cdef class PyDiffusionProblem:
    cdef DiffusionProblem *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, num_elements):
        self.thisptr = new DiffusionProblem(num_elements)
    def __dealloc__(self):
        del self.thisptr
    def read_dirichlet_indices(self, indices):
        self.thisptr.read_dirichlet_indices(indices)
    def return_zero_vector(self):
        return self.thisptr.return_zero_vector()
    def matrix_vector_multiply(self, vec):
        return self.thisptr.matrix_vector_multiply(vec)
    def size_of_system(self):
        return self.thisptr.size_of_system()
    def plot_solution(self, vec):
        self.thisptr.plot_solution(vec)
        
