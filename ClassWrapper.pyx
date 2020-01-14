# distutils: language=c++

cdef class PyDiffusionProblem:
  cdef DiffusionProblem *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self, num_elements):
    self.thisptr = new DiffusionProblem(num_elements, num_elements, 1.)
  def __dealloc__(self):
    del self.thisptr
  def read_dirichlet_indices(self, indices):
    self.thisptr.read_dirichlet_indices(indices)
  def return_zero_vector(self):
    return self.thisptr.return_zero_vector()
  def matrix_vector_multiply(self, vec):
    return self.thisptr.matrix_vector_multiply(vec)
  def size_of_system(self):
    return self.thisptr.size_of_system()
  def plot_option(self, option, value):
    if isinstance(option,str): # Python3 version - use unicode for Python 2
      option = option.encode()
    if isinstance(value,str): # Python3 version - use unicode for Python 2
      value = value.encode()
    return_val = self.thisptr.plot_option(option, value)
    if isinstance(return_val,bytes): # Python3 version - use unicode for Python 2
      return_val = return_val.decode()
    return return_val
  def plot_solution(self, vec):
    self.thisptr.plot_solution(vec)

cdef class PyElasticityProblem:
  cdef ElasticityProblem *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self, num_elements):
    self.thisptr = new ElasticityProblem(num_elements, num_elements, 1.)
  def __dealloc__(self):
    del self.thisptr
  def read_dirichlet_indices(self, indices):
    self.thisptr.read_dirichlet_indices(indices)
  def return_zero_vector(self):
    return self.thisptr.return_zero_vector()
  def matrix_vector_multiply(self, vec):
    return self.thisptr.matrix_vector_multiply(vec)
  def size_of_system(self):
    return self.thisptr.size_of_system()
  def plot_option(self, option, value):
    if isinstance(option,str): # Python3 version - use unicode for Python 2
      option = option.encode()
    if isinstance(value,str): # Python3 version - use unicode for Python 2
      value = value.encode()
    return_val = self.thisptr.plot_option(option, value)
    if isinstance(return_val,bytes): # Python3 version - use unicode for Python 2
      return_val = return_val.decode()
    return return_val
  def plot_solution(self, vec):
    self.thisptr.plot_solution(vec)
