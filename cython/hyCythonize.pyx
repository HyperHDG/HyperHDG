# distutils: language=c++

def hyPyCythonize(names, py_ver_maj, py_ver_min):
  if isinstance(names[0],str):
    names = [i.encode() for i in names]
  filename = hyCythonize(names, py_ver_maj, py_ver_min)
  if isinstance(filename,bytes):
    filename = filename.decode()
  return filename
