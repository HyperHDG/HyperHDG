import os
import importlib

def importer(names):
  if not os.path.isfile("build/cythonize"):
    os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from hyImporter import hyImport
  retval = hyImport(names)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
