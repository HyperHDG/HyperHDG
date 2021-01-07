## Transform C++ class name to cython file name.
def cython_from_cpp(name):
  helper = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
  return re.sub('([a-z0-9])([A-Z])', r'\1_\2', helper).lower()

## \brief   Generate names of auxiliary classes that will be constructed (internal use only).
def files(conf):
  cpp_class    = "GlobalLoop::" + conf.global_loop + "<Topology::" + conf.topology \
                 + ",Geometry::" + conf.geometry + ",NodeDescriptor::" \
                 + conf.node_descriptor + ",LocalSolver::" + conf.local_solver + " >"
  cython_file  = cpp_class + "_deb" if conf.debug_mode else cpp_class + "_rel"
  cython_file  = re.sub(' ', '', cython_file)
  cython_file  = re.sub('\+|\-|\*|\/|<|>|\,|\:', '_', cython_file)
  cython_class = cython_file + "CP"
  return cpp_class, cython_file, cython_class