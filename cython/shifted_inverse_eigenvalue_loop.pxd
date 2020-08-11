# C++: HyperHDG/global_loop/shifted_inverse_eigenvalue.hxx

# CyReplace_Number: 4      // For setting default values for CyReplace, their amount must be known.
# CyReplace03: double      // Set default value for CyReplace03 to "double".
# CyReplace04: double      // Set default value for CyReplace04 to "double".

from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/shifted_inverse_eigenvalue.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    vector[unsigned int] dirichlet_indices ()
    vector[ CyReplace03 ] matrix_vector_multiply (vector[ CyReplace03 ], CyReplace03 )
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
