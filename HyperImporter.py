import os
from importlib import import_module

def importer(names):
  if not os.path.isfile("build/cythonize"):
    os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from importer import importing
  retval = importing(names)
  mod = import_module(retval)
  return getattr(mod, retval)
