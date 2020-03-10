import os
from os import path

def importer(names):
#   if not path.isfile("build/cythonize"):
  os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize; build/cythonize")
  from importer import importing
  importing(names)
  from PyDiffusionProblem import PyDiffusionProblem
  return PyDiffusionProblem
  print("Importer is running!")
