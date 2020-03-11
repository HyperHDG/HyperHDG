import os
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
#  build directory is not present, the Python call is invalid and the program will construct both
#  needed directories and terminate, leaving the operating system in a consistent state for the
#  Python code to be restarted. If this has beed done, the program compiles the Cythonize.cxx
#  file that does all the work: It takes a .pyx and a .pxd file and formats them to allow the
#  coupling of Python and C++ code and does also all compilation steps. The file can be started from
#  Python using the hyCythonizer module (which is generated by compilation of Cythonize.cxx).
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
  if not os.path.isfile("build/cythonize"):
    if not os.path.exists("build/CythonFiles"):
      print("\nFull build path is not available. The program corrects this and eventually dies.")
      print("Afterwards, it can be restarted and should work out!\n")
    os.system("mkdir -p build build/CythonFiles build/SharedObjects")
    os.system("g++ Cythonize.cxx -std=c++17 -Iinclude -o build/cythonize -lstdc++fs")
    os.system("./build/cythonize")
  from hyCythonizer import hyPyCythonize
  retval = hyPyCythonize(names)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
