cdef class PyDiffusionProblem:
  cdef DiffusionProblem *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self, num_elements):
    self.thisptr = new DiffusionProblem(num_elements)
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
  def plot_solution(self, vec):
    self.thisptr.plot_solution(vec)
