# C++: HyperHDG/AbstractProblem.hxx    (This line allows for just in time compile optimization!)

# CyReplace_Number: 3      // For setting default values for CyReplace, their amount must be known.
# CyReplace03: double      // Set default value for CyReplace03 to "double".

from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/AbstractProblem.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , double) except +
    CythonClassName ( CyReplace01 , double) except +
    CythonClassName ( CyReplace01 ) except +
    void read_dirichlet_indices (vector[unsigned int])
    vector[ CyReplace03 ] return_zero_vector ()
    vector[ CyReplace03 ] matrix_vector_multiply (vector[ CyReplace03 ])
    vector[ CyReplace03 ] assemble_rhs (vector[ CyReplace03 ])
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ])
