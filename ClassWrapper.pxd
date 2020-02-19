from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "AbstractProblem.hxx":
  cdef cppclass DiffusionProblemNaive "DiffusionProblemRegularNaive<2,3,1>":
    DiffusionProblemNaive(vector[unsigned int], vector[unsigned int], double) except +
    void read_dirichlet_indices(vector[unsigned int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])

cdef extern from "AbstractProblem.hxx":
  cdef cppclass DiffusionProblemNaiveF "DiffusionProblemRegularNaiveF<2,3,1>":
    DiffusionProblemNaiveF(vector[unsigned int], vector[unsigned int], float) except +
    void read_dirichlet_indices(vector[unsigned int])
    vector[float] return_zero_vector[float]()
    vector[float] matrix_vector_multiply(vector[float])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[float])

cdef extern from "AbstractProblem.hxx":
  cdef cppclass ElasticityProblem "ElasticityProblemFile<1,2,1>":
    ElasticityProblem(string, double) except +
    void read_dirichlet_indices(vector[unsigned int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])
"""
cdef extern from "AbstractProblem.hxx":
  cdef cppclass ElasticityBBeam "ElasticityBBeam<1,2,1>":
    ElasticityBBeam(string, double) except +
    void read_dirichlet_indices(vector[unsigned int])
    vector[double] return_zero_vector()
    vector[double] matrix_vector_multiply(vector[double])
    int size_of_system()
    string plot_option(string, string)
    void plot_solution(vector[double])
"""
