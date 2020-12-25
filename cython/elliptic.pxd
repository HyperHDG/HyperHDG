from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/elliptic.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    void read_dirichlet_indices (vector[unsigned int])
    vector[ CyReplace03 ] return_zero_vector ()
    vector[ CyReplace03 ] matrix_vector_multiply (vector[ CyReplace03 ])
    vector[ CyReplace03 ] total_flux_vector (vector[ CyReplace03 ])
    CyReplace03 calculate_L2_error (vector[ CyReplace03 ])
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ])