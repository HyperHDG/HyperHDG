import os
import importlib

def hyImport(names):
  if not os.path.isfile("build/cythonize"):
    if not os.path.exists("build/CythonFiles"):
      print("\nFull build path is not available. The program corrects this and eventually dies.")
      print("Afterwards, it can be restarted and should work out!\n")
    os.system("mkdir -p build build/CythonFiles")
    os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from hyCythonizer import hyPyCythonize
  retval = hyPyCythonize(names)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
