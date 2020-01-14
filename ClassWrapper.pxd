from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "AbstractProblem.h":
  cdef cppclass DiffusionProblem "DiffusionProblemRegular<1,3,1>":
    DiffusionProblem(vector[int], vector[int], double) except +
    void read_dirichlet_indices(vector[int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])

cdef extern from "AbstractProblem.h":
  cdef cppclass ElasticityProblem "ElasticityProblemRegular<1,3,1>":
    ElasticityProblem(vector[int], vector[int], double) except +
    void read_dirichlet_indices(vector[int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])
