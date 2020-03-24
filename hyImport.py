import os
import sys
import importlib

## \brief   Function to import classes of the HyperHDG package using Cython.
#
#  This function takes a vector / list of at least two strings. The first string refers to the name
#  of the .pyx and .pxd files located in the cython directory. Each of the .pxd files is clearly
#  associated to one and only one C++ file (.hxx or .cxx). The second string refers to the class
#  that is to be imported (and is defined in the C++ file). Further strings might refer to embedded
#  template parameters. This has, however, not yet been fully implemented.
#
#  First, the method checks whether the build and build/CythonFiles directories are present. If the
#  build directory is not present, the Python program will construct both  needed directories. If 
#  this has beed done, the program compiles the Cythonize.cxx file that does all the work:
#  It takes  a .pyx and a .pxd file and formats them to allow the coupling of Python and C++ code
#  and does also all compilation steps. The file can be started from Python using the hyCythonizer
#  module (which is generated by compilation of Cythonize.cxx).
#
#  Afterwards, the program imports the class indicated by the vector of names and returns this class
#  to the Python script that started the function. By this, the class is also imported to that
#  script.
#
#  \param   names   Vector containing names of class to be imported.
#  \retval  class   The C++ class (not class member) to be used in the Python code.
#
#  \authors   Guido Kanschat, University of Heidelberg, 2020.
#  \authors   Andreas Rupp, University of Heidelberg, 2020.
def hyImport(names):
  ver_major = sys.version_info.major
  ver_minor = sys.version_info.minor

  if not os.path.isfile("build/cythonize"):
    if not os.path.isfile("/usr/include/python"+str(ver_major)+"."+str(ver_minor)+"/Python.h"):
      print("\nThe current Python version seems not to have an include file.\n")
      print("This will most likely result in an error!\n")
      print("Check your Python version to be not of m or dm type which is not fully supported.\n\n")
    os.system("mkdir -p build build/CythonFiles build/SharedObjects")
    os.system("g++ Cythonize.cxx -DPYVERMAJ=" + str(ver_major) + " -DPYVERMIN=" + str(ver_minor) +
              " -std=c++17 -Iinclude -o build/cythonize -lstdc++fs")
    os.system("./build/cythonize")

  try:
    from hyCythonizer import hyPyCythonize
  except ImportError as error:
    sys.path.append(os.path.dirname(__file__) + "/build/SharedObjects")
    from hyCythonizer import hyPyCythonize

  retval = hyPyCythonize(names, ver_major, ver_minor)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
