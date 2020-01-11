from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "DiffusionProblem.h":
  cdef cppclass DiffusionProblem "DiffusionProblemRegular<1,3,1>":
    DiffusionProblem(vector[int], vector[int]) except +
    void read_dirichlet_indices(vector[int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])
