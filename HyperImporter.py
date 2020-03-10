import os
from os import path
import importlib

def importer(names):
#   if not path.isfile("build/cythonize"):
  os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from importer import importing
  retval = importing(names)
  from DiffusionProblemRegularNaive_1_3_1_ import DiffusionProblemRegularNaive_1_3_1_
  return DiffusionProblemRegularNaive_1_3_1_
  print("Importer is running!")
