# distutils: language=c++

cdef class PythonClassName :
  cdef CythonClassName *thisptr # hold a C++ instance which we're wrapping
  def __cinit__(self, topo_constr, geom_constr = 'default', lsol_constr = 'default'):
    if geom_constr == 'default':
      if isinstance(topo_constr,str):
        topo_constr = topo_constr.encode()
      if lsol_constr == 'default':
        self.thisptr = new CythonClassName (topo_constr)
      else:
        self.thisptr = new CythonClassName (topo_constr, lsol_constr)
    else:
      if isinstance(topo_constr,str):
        topo_constr = topo_constr.encode()
      if isinstance(geom_constr,str):
        geom_constr = geom_constr.encode()
      if lsol_constr == 'default':
        print("Python defaults lsol_constr to 1 and might overwrite C++ settings at this point.\n")
        # self.thisptr = new CythonClassName (topo_constr, geom_constr, 1.)
      else:
        self.thisptr = new CythonClassName (topo_constr, geom_constr, lsol_constr)
  def __dealloc__(self):
    del self.thisptr
  def dirichlet_indices(self):
    return self.thisptr.dirichlet_indices ()
  def trace_to_flux(self, vec, time = 0.):
    return self.thisptr.trace_to_flux (vec, time)
  def size_of_system(self):
    return self.thisptr.size_of_system ()
  def plot_option(self, option, value):
    if isinstance(option,str):
      option = option.encode()
    if isinstance(value,str):
      value = value.encode()
    return_val = self.thisptr.plot_option (option, value)
    if isinstance(return_val,bytes):
      return_val = return_val.decode()
    return return_val
  def plot_solution(self, vec, time = 0.):
    self.thisptr.plot_solution (vec, time)
  # Condensed matrix in coordinate format (entries with equal indices have to be summed up).
  # factorize_once = True computes the local matrix only once per combination of boundary
  # conditions and reuses it for all hyperedges. This is only correct for uniform meshes with
  # coefficients that are constant in space (not, e.g., for space dependent diffusion tensors).
  # Violations are detected for a sample of hyperedges and raise a ValueError.
  def sparse_stiff_mat(self, sigma = 0., factorize_once = False):
    helper = self.thisptr.trace_to_flux_mat(sigma, factorize_once)
    return helper.get_cols(), helper.get_rows(), helper.get_values()
  # Let trace_to_flux use a sparse matrix that is assembled once (faster if trace_to_flux is called
  # often, needs memory for the matrix). The matrix is reassembled if trace_to_flux is called with
  # another shift, unless shift_independent = True.
  def cache_operator(self, enable = True, factorize_once = False, shift_independent = False):
    self.thisptr.cache_operator(enable, factorize_once, shift_independent)
  def refine(self, n_ref = 'default'):
    if n_ref == 'default':
      return self.thisptr.get_refinement()
    else:
      self.thisptr.set_refinement(n_ref)
      return n_ref
