from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/mass_approx_eigenvalue.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    vector[unsigned int] dirichlet_indices ()
    vector[ CyReplace03 ] trace_to_flux (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] trace_to_mass_flux (vector[ CyReplace03 ], CyReplace03 )
    int size_of_system ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
    unsigned int get_refinement()
    void set_refinement(unsigned int)
