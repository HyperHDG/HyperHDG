# distutils: language=c++

def hyPyCythonize(names, filenames, py_ver_maj, py_ver_min, force_comp):
  if isinstance(names[0],str):
    names = [i.encode() for i in names]
  if isinstance(filenames[0],str):
    filenames = [i.encode() for i in filenames]
  classname = hyCythonize(names, filenames, py_ver_maj, py_ver_min, force_comp)
  if isinstance(classname,bytes):
    classname = classname.decode()
  return classname
