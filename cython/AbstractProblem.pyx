# distutils: language=c++

cdef class PythonClassName :
  cdef CythonClassName *thisptr # hold a C++ instance which we're wrapping
  def __cinit__(self, topo_constr, geom_constr = 'default', lsol_constr = 'default'):
    if geom_constr == 'default':
      if isinstance(topo_constr,str): # Python3 version - use unicode for Python 2
        topo_constr = topo_constr.encode()
      if lsol_constr == 'default':
        self.thisptr = new CythonClassName (topo_constr)
      else:
        self.thisptr = new CythonClassName (topo_constr, lsol_constr)
    else:
      if isinstance(topo_constr,str): # Python3 version - use unicode for Python 2
        topo_constr = topo_constr.encode()
      if isinstance(geom_constr,str): # Python3 version - use unicode for Python 2
        geom_constr = geom_constr.encode()
      if lsol_constr == 'default':
        print("Python defaults lsol_constr to 1 and might overwrite C++ settings at this point.\n")
        #self.thisptr = new CythonClassName (topo_constr, geom_constr, 1.)
      else:
        self.thisptr = new CythonClassName (topo_constr, geom_constr, lsol_constr)
  def __dealloc__(self):
    del self.thisptr
  def read_dirichlet_indices(self, indices):
    self.thisptr.read_dirichlet_indices (indices)
  def dirichlet_indices(self):
    return self.thisptr.dirichlet_indices ()
  def return_zero_vector(self):
    return self.thisptr.return_zero_vector[ CyReplace03 ] ()
  def matrix_vector_multiply(self, vec, time = 0.):
    return self.thisptr.matrix_vector_multiply (vec, time)
  def matrix_vector_der_multiply(self, vec, eig, vec_val, eig_val):
    return self.thisptr.matrix_vector_der_multiply (vec, eig, vec_val, eig_val)
  def total_flux_vector(self, vec, time = 0.):
    return self.thisptr.total_flux_vector (vec, time)
  def set_data(self, vec, time = 0.):
    self.thisptr.set_data (vec, time)
  def initial_flux_vector(self, vec, time = 0.):
    return self.thisptr.initial_flux_vector (vec, time)
  def mass_matrix_multiply(self, vec, time = 0.):
    return self.thisptr.mass_matrix_multiply (vec, time)
  def total_mass_vector(self, vec, time = 0.):
    return self.thisptr.total_mass_vector (vec, time)
  def calculate_L2_error(self, vec, time = 0.):
    return self.thisptr.calculate_L2_error (vec, time)
  def calculate_L2_error_temp(self, vec_new, vec_old, delta_time, time):
    return self.thisptr.calculate_L2_error_temp (vec_new, vec_old, delta_time, time)
  def size_of_system(self):
    return self.thisptr.size_of_system ()
  def plot_option(self, option, value):
    if isinstance(option,str): # Python3 version - use unicode for Python 2
      option = option.encode()
    if isinstance(value,str): # Python3 version - use unicode for Python 2
      value = value.encode()
    return_val = self.thisptr.plot_option (option, value)
    if isinstance(return_val,bytes): # Python3 version - use unicode for Python 2
      return_val = return_val.decode()
    return return_val
  def plot_solution(self, vec, time = 0.):
    self.thisptr.plot_solution (vec, time)
