# distutils: language=c++

def importing(names):
  if isinstance(names[0],str):
    names = [i.encode() for i in names]
  name = hyCython(names)
  if isinstance(name,bytes):
    name = name.decode()
  return name
