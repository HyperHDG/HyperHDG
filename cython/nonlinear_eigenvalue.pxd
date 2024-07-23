from libcpp.vector cimport vector
from libcpp.string cimport string

IncludeFiles

cdef extern from "<HyperHDG/global_loop/prototype.hxx>" :
  cdef cppclass sparse_mat "sparse_mat<std::vector< CyReplace03 > >" :
    sparse_mat ( unsigned int ) except +
    sparse_mat () except +
    vector[ unsigned int ] get_cols()
    vector[ unsigned int ] get_rows()
    vector[ CyReplace03 ] get_values()

cdef extern from "<HyperHDG/global_loop/nonlinear_eigenvalue.hxx>" :
  cdef cppclass CythonClassName C++ClassName :
    CythonClassName ( CyReplace01 , CyReplace02 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 , CyReplace04 ) except +
    CythonClassName ( CyReplace01 ) except +
    vector[ CyReplace03 ] trace_to_flux (vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] jacobian_of_trace_to_flux (vector[ CyReplace03 ], CyReplace03 , vector[ CyReplace03 ], CyReplace03 )
    vector[ CyReplace03 ] make_initial (vector[ CyReplace03 ], CyReplace03 )
    int size_of_system ()
    vector[ unsigned int ] dirichlet_nodes ()
    string plot_option (string, string)
    void plot_solution (vector[ CyReplace03 ], CyReplace03 )
    unsigned int get_refinement ()
    void set_refinement (unsigned int)
    sparse_mat trace_to_flux_mat( CyReplace03 )
    sparse_mat jacobian_of_trace_to_flux_mat( vector[ CyReplace03 ], CyReplace03 )
