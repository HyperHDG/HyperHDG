# C++: HyperHDG/AbstractProblem.hxx    (This line allows for just in time compile optimization!)

# CyReplace_Number: 4      // For setting default values for CyReplace, their amount must be known.
# CyReplace03: double      // Set default value for CyReplace03 to "double".
# CyReplace04: double      // Set default value for CyReplace04 to "double".

from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/AbstractProblem.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    void read_dirichlet_indices (vector[unsigned int])
    vector[ CyReplace03 ] return_zero_vector[ CyReplace03 ] ()
    vector[ CyReplace03 ] matrix_vector_multiply (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] total_flux_vector (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] initial_flux_vector (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] mass_matrix_multiply (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] total_mass_vector (vector[ CyReplace03 ], CyReplace03 )
    CyReplace03 calculate_L2_error (vector[ CyReplace03 ], CyReplace03 )
    CyReplace03 calculate_L2_error_temp (vector[ CyReplace03 ], vector[ CyReplace03 ], CyReplace03 , CyReplace03 )
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
