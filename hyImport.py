import os
import importlib

def hyImport(names):
  if not os.path.isfile("build/cythonize"):
    os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from hyCythonizer import hyPyCythonize
  retval = hyPyCythonize(names)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
