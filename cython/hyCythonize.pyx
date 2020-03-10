# distutils: language=c++

def hyImport(names):
  if isinstance(names[0],str):
    names = [i.encode() for i in names]
  filename = hyCythonize(names)
  if isinstance(filename,bytes):
    filename = filename.decode()
  return filename
