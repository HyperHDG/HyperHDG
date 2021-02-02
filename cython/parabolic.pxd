from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/parabolic.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    vector[ CyReplace03 ] zero_vector ()
    vector[ CyReplace03 ] trace_to_flux (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] residual_flux (vector[ CyReplace03 ], CyReplace03 )
    void set_data (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] make_initial (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] errors (vector[ CyReplace03 ], CyReplace03 )
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
    unsigned int get_refinement()
    void set_refinement(unsigned int)
