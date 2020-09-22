import os, sys, importlib

## \brief   Function to import classes of the HyperHDG package using Cython.
#
#  This function takes a vector / list of at least two strings. The first string refers to the name
#  of the .pyx and .pxd files located in the cython directory. Each of the .pxd files is clearly
#  associated to one and only one C++ file (.hxx or .cxx). The second string refers to the class
#  that is to be imported (and is defined in the C++ file). Further strings might refer to embedded
#  template parameters or other strings that need to be substituted. Thus, if a .pxd file contains
#  the key word \c CyReplace followed by a two digit number (witout space between word and number),
#  it is replace by the \c vector name's number + 2 entry (the first two names have predefined
#  meaning, cf. above).
#
#  Additionally, a vector of filenames can be specified. These filenames refer to files that need
#  to be imported/included to allow the C++ classes to be compilable. For example, if the class
#  AbstractProblem is to be cythonized it receives a Topology, a Geometry, and a Local Solver as
#  template parameters. The files in which the specified Topology, Geometry, and Local Solver are
#  defined need to be available to the compiler and therefore be added in the filelist.
#
#  First, the method checks whether the build and build/cython_files directories are present. If the
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
#  \param   names       Vector containing names of class to be imported.
#  \param   filenames   Vector containing names of files that need to be included. Default is empty.
#  \param   force_comp  Force recompilation of the C++ class. Default is False.
#  \param   debug_mode  Compile library in debub mode (not in release mode). Default is False.
#  \retval  class       The C++ class (not class member) to be used in the Python code.
#
#  \authors   Guido Kanschat, Heidelberg University, 2020.
#  \authors   Andreas Rupp, Heidelberg University, 2020.
def cython_import(names, filenames = [], force_comp = False, debug_mode = False):
  ver_major = sys.version_info.major
  ver_minor = sys.version_info.minor

  if not os.path.isfile("build/cythonize"):
    if not os.path.isfile("/usr/include/python"+str(ver_major)+"."+str(ver_minor)+"/Python.h"):
      print("\nThe current Python version seems not to have an include file.\n")
      print("This will most likely result in an error!\n")
      print("Check your Python version to be not of m or dm type which is not fully supported.\n\n")
    os.system("mkdir -p build build/cython_files build/shared_objects output")
    os.system("g++ cython/cythonize.cxx -DPYVERMAJ=" + str(ver_major) + " -DPYVERMIN=" +
              str(ver_minor) + " -std=c++17 -Iinclude -o build/cythonize -lstdc++fs")
    os.system("./build/cythonize")

  try:
    from hyper_cythonizer import hyper_cythonize
  except ImportError as error:
    sys.path.append(os.path.dirname(__file__) + "/build/shared_objects")
    from hyper_cythonizer import hyper_cythonize

  retval = hyper_cythonize(names, filenames, ver_major, ver_minor, force_comp, debug_mode)
  mod = importlib.import_module(retval)
  return getattr(mod, retval)
