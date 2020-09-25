# distutils: language=c++

def hyper_cythonize(names, filenames, py_ver_maj, py_ver_min, debug_mode):
  if isinstance(names[0],str):
    names = [i.encode() for i in names]
  if isinstance(filenames[0],str):
    filenames = [i.encode() for i in filenames]
  classname = cythonize(names, filenames, py_ver_maj, py_ver_min, debug_mode)
  if isinstance(classname,bytes):
    classname = classname.decode()
  return classname
