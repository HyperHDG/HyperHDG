# C++: AbstractProblem.hxx    (This line allows for just in time compile optimization!)

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<Geom_UnitCube.hxx>" : pass
cdef extern from "<HyperHDG/LocalSolver/Diffusion.hxx>" : pass

cdef extern from "<AbstractProblem.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName (vector[unsigned int], vector[unsigned int], double) except +
    void read_dirichlet_indices (vector[unsigned int])
    vector[double] return_zero_vector ()
    vector[double] matrix_vector_multiply (vector[double])
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[double])
