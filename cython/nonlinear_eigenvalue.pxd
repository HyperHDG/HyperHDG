from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/nonlinear_eigenvalue.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    vector[ CyReplace03 ] trace_to_flux (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] jacobian_of_trace_to_flux (vector[ CyReplace03 ], CyReplace03 , vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] make_initial (vector[ CyReplace03 ], CyReplace03 )
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
