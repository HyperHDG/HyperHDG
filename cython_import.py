## \package cython_import
# 
#  \brief   Python interface for HyperHDG library.
#
#  This file is the main ingredient to use the HyperHDG library within Python scripts. Its function
#  \c cython_import compiles the necessary components just-in-time and provides the as classes that
#  can be used within Python. To this end, the Cython package is heavily used. 
#
#  \authors   Guido Kanschat, Heidelberg University, 2020.
#  \authors   Andreas Rupp, Heidelberg University, 2020.

## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
COMPILE_COM = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
COMPILE_INC = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
COMPILE_FLG = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
COMPILE_STD = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
LINK_COM = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
LINK_FLG = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
LINK_LIB = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
CYTHON_COM = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
CYTHON_FLG = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
PY_VER_MAJ = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
PY_VER_MIN = ""
## \brief   Parameter that might be set by cmake using the function get_cmakes() (cf. below).
PYTHON_DIR = ""

import configparser, os, sys, importlib, glob, re, datetime

## \brief   Object that comprises all information for HyperHDG to create a problem.
#
#  \authors   Guido Kanschat, Heidelberg University, 2020.
#  \authors   Andreas Rupp, Heidelberg University, 2020.
class hyperhdg_constructor:
  global_loop         = ""
  local_solver        = ""
  topology            = ""
  geometry            = ""
  node_descriptor     = ""
  cython_replacements = []
  include_files       = []
  debug_mode          = False
  allow_file_output   = True
  def is_consistent(self):
    if not (isinstance(self.global_loop, str) and self.global_loop != ""):
      return False
    if not (isinstance(self.local_solver, str) and  self.local_solver != ""):
      return False
    if not (isinstance(self.topology, str) and self.topology != ""):
      return False
    if not (isinstance(self.geometry, str) and self.geometry != ""):
      return False
    if not (isinstance(self.node_descriptor, str) and self.node_descriptor != ""):
      return False
    if not isinstance(self.debug_mode, bool):
      return False
    if not isinstance(self.allow_file_output, bool):
      return False
    if not isinstance(self.cython_replacements, list):
      return False
    for entry in self.cython_replacements:
      if not isinstance(entry, str):
        return False
    if not isinstance(self.include_files, list):
      return False
    for entry in self.include_files:
      if not isinstance(entry, str):
        return False
    return True

## \brief   Function to import classes of the HyperHDG package using Cython.
#
#  \param   constructor hyperhdg_constructor object.
#
#  \authors   Guido Kanschat, Heidelberg University, 2020.
#  \authors   Andreas Rupp, Heidelberg University, 2020.
def cython_import(constructor):
  start_time = datetime.datetime.now()
  # Check that constructor is appropriatelly filled.
  assert isinstance(constructor, hyperhdg_constructor) and constructor.is_consistent()

  # Start program.
  get_cmakes()
  print("Cythonizing ... ", end='', flush=True)

  # Create folders and log files and check for consistency.
  os.system("mkdir -p " + main_path() + "/build " + main_path() + "/build/cython_files " \
            + main_path() + "/build/shared_objects")
  if not os.path.isfile(main_path() + "/build/cython_log.txt"):
    file = open(main_path() + "/build/cython_log.txt", "w")
    file.write("Python version: " + str(sys.version_info))
    file.close()
  else:
    file = open(main_path() + "/build/cython_log.txt", "r")
    assert file.readline() == "Python version: " + str(sys.version_info), \
           "Wrong Python version!\n \
            Remove build directory and reinstall HyperHDG or use previous Python version."
    file.close()

  # Evaluate configfile, define include files, and define dependent files.
  cy_replace      = evaluate_config(constructor)
  include_string  = extract_includes(constructor)
  cpp_class, python_class, cython_class = file_names(constructor)

  compilation_necessary = need_compile(constructor, python_class)
  if compilation_necessary:
    # Copy pyx and pxd files from cython directory.
    for file_end in ["pxd", "pyx"]:
      with open(main_path() + "/cython/" + cython_from_cpp(constructor.global_loop) + "." \
        + file_end, "r") as file:
        content = file.read()
      content = re.sub("C\+\+ClassName", "\"" + cpp_class + "\"", content)
      content = re.sub("CythonClassName", cython_class, content)
      content = re.sub("PythonClassName", python_class, content)
      content = re.sub("IncludeFiles", include_string, content)
      for i in range(len(cy_replace)):
        content = re.sub("CyReplace" + '%02d' % (i+1), cy_replace[i], content)
      with open(main_path() + "/build/cython_files/" + python_class + "." + file_end, "w") as file:
        file.write(content)
    # Prepare the compilation commands.
    cython_command, compile_command, link_command = get_commands(python_class)
    if not (constructor.debug_mode):
      compile_command += " -DNDEBUG";
    if not (constructor.allow_file_output):
      compile_command += " -DNOFILEOUT"
    #Actually compile the prepared files.
    assert os.system(cython_command) == 0
    assert os.system(compile_command) == 0
    assert os.system(link_command) == 0

  try:
    mod = importlib.import_module(python_class)
  except ImportError as error:
    sys.path.append(main_path() + "/build/shared_objects")
    mod = importlib.import_module(python_class)

  delta_time_ms = 1000 * (datetime.datetime.now() - start_time).total_seconds()
  if compilation_necessary:
    print("DONE with compilation in " + "{:.2f}".format(delta_time_ms) + " milliseconds.")
  else:
    print("DONE without compilation in " + "{:.2f}".format(delta_time_ms) + " milliseconds.")

  return getattr(mod, python_class)

## \brief   Extract classname from name that might contain template arguemnts.
def extract_classname(fullname):
  index = re.search('\<|\ ', fullname)
  return fullname[0:index.start() if index != None else len(fullname)]

## \brief   Find file with definition of classname in directory folder.
def find_definition(folder, classname):
  for file in glob.glob(main_path() + "/include/HyperHDG/" + folder + "/*.hxx"):
    with open(file, "r") as hxxfile:
      if ("class " + classname or "struct " + classname) in hxxfile.read():
        return re.sub(main_path() + "/include/", '', file)
  assert False, "File containing defintion of " + classname + " has not been found!"

## \brief   Evaluate config file that needs to be present for all .pyx/.pxd files.
def evaluate_config(constructor):
  config = configparser.ConfigParser()
  config.read(main_path() + "/cython/" + cython_from_cpp(constructor.global_loop) + ".cfg")
  n_replacements = int(config['default']['n_replacements'])
  cy_replace = []
  for i in range(n_replacements):
    if i < len(constructor.cython_replacements) and constructor.cython_replacements[i] != "":
      cy_replace.append(constructor.cython_replacements[i])
    else:
      assert config.has_option('default','replacement' + '%02d' % (i+1))
      cy_replace.append(config['default']['replacement' + '%02d' % (i+1)])
  assert n_replacements == len(cy_replace)
  return cy_replace

## \brief   Read out parameters for Cython compilation provided by CMAKE.
def get_cmakes():
  global COMPILE_COM, COMPILE_INC, COMPILE_FLG, COMPILE_STD, LINK_COM, LINK_FLG, LINK_LIB, \
    CYTHON_COM, CYTHON_FLG, PY_VER_MAJ, PY_VER_MIN, PYTHON_DIR
  if not os.path.isfile(main_path() + "/build/cmake_cython.cfg"):
    print("CMAKE files do not exist ... using default values for Cython!")
    COMPILE_COM = "g++"
    COMPILE_INC = "-I. -Iinclude -Isubmodules/tensor_product_chain_complex.git/include \
      -I/usr/include/python" + str(sys.version_info.major) + "." + str(sys.version_info.minor)
    COMPILE_FLG = "-pthread -g -fwrapv -O2 -Wall -pedantic -g -fstack-protector-strong -Wformat \
      -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC"
    COMPILE_STD = "17"
    LINK_COM = "g++"
    LINK_FLG = "-pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
      -Wl,-z,relro -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong -Wformat \
      -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2"
    LINK_LIB = "-llapack"
    CYTHON_COM = "cython3"
    CYTHON_FLG = "-3 --cplus"
    PY_VER_MAJ = str(sys.version_info.major)
    PY_VER_MIN = str(sys.version_info.minor)
    PYTHON_DIR = "/usr/include/python"+str(sys.version_info.major)+"."+str(sys.version_info.minor)
  else:
    config = configparser.ConfigParser()
    config.read(main_path() + "/build/cmake_cython.cfg")
    COMPILE_COM = config['compiler']['command']
    COMPILE_INC = config['compiler']['includes']
    COMPILE_FLG = config['compiler']['flags']
    COMPILE_STD = config['compiler']['standard']
    LINK_COM = config['linker']['command']
    LINK_FLG = config['linker']['flags']
    LINK_LIB = config['linker']['libs']
    CYTHON_COM = config['cython']['command']
    CYTHON_FLG = config['cython']['flags']
    PY_VER_MAJ = config['python']['ver_maj']
    PY_VER_MIN = config['python']['ver_min']
    PYTHON_DIR = config['python']['dir']
  if int(PY_VER_MAJ) < 3:
    print("Python versions below 3 are not supported!")
  assert PYTHON_DIR in COMPILE_INC, "Python directory must be included!"
  assert int(PY_VER_MAJ) == sys.version_info.major and int(PY_VER_MIN) == sys.version_info.minor, \
         "Utilized Python version is not CMAKE's Python version!"

## \brief   Add the files that need to be included for defining the problem to include list.
def extract_includes(constructor):
  helper = find_definition("geometry", extract_classname(constructor.geometry))
  if helper not in constructor.include_files:
    constructor.include_files.append(helper)
  helper = find_definition("node_descriptor", extract_classname(constructor.node_descriptor))
  if helper not in constructor.include_files:
    constructor.include_files.append(helper)
  helper = find_definition("topology", extract_classname(constructor.topology))
  if helper not in constructor.include_files:
    constructor.include_files.append(helper)
  helper = find_definition("global_loop", extract_classname(constructor.global_loop))
  if helper not in constructor.include_files:
    constructor.include_files.append(helper)
  helper = find_definition("local_solver", extract_classname(constructor.local_solver))
  if helper not in constructor.include_files:
    constructor.include_files.append(helper)
  include_string = ""
  for x in constructor.include_files:
    include_string = include_string + "cdef extern from \"<" + x + ">\": pass\n"
  return include_string

## \brief   Generate names of auxiliary classes that will be constructed (internal use only).
def file_names(constructor):
  cpp_class    = "GlobalLoop::" + constructor.global_loop + "<Topology::" + constructor.topology \
                 + ",Geometry::" + constructor.geometry + ",NodeDescriptor::" \
                 + constructor.node_descriptor + ",LocalSolver::" + constructor.local_solver + " >"
  cython_file  = cpp_class + "_deb" if constructor.debug_mode else cpp_class + "_rel"
  cython_file  = re.sub(' ', '', cython_file)
  cython_file  = re.sub('\+|\-|\*|\/|<|>|\,|\:', '_', cython_file)
  cython_class = cython_file + "CP"
  return cpp_class, cython_file, cython_class

## \brief   Generate shell commands from CMAKE parameters.
def get_commands(python_class):
  global COMPILE_COM, COMPILE_INC, COMPILE_FLG, COMPILE_STD, LINK_COM, LINK_FLG, LINK_LIB, \
    CYTHON_COM, CYTHON_FLG
  cython_command = "cd " + main_path() + "/build/cython_files/; " \
    + CYTHON_COM + " " + CYTHON_FLG + " " + python_class + ".pyx";
  compile_command = "cd " + main_path() + "; " \
    + COMPILE_COM + " " + COMPILE_INC + " " + COMPILE_FLG + " --std=gnu++" + COMPILE_STD + " -c \
    ./build/cython_files/" + python_class + ".cpp -o ./build/cython_files/" + python_class + ".o";
  link_command = "cd " + main_path() + "; " \
    + LINK_COM + " " + LINK_FLG + " ./build/cython_files/" + python_class + ".o " \
    + " -o build/shared_objects/" + python_class + ".so " + LINK_LIB;
  return cython_command, compile_command, link_command

## \brief   Check whether recompilation of executable is necessary.
def need_compile(constructor, python_class):
  global PYTHON_DIR
  if not os.path.isfile(main_path() + "/build/shared_objects/" + python_class + ".so"):
    return True
  time_so = os.stat(main_path() + "/build/shared_objects/" + python_class + ".so").st_mtime
  assert os.path.isfile(PYTHON_DIR + "/Python.h"), "Include file for Python needs to exist!"
  if time_so < os.stat(PYTHON_DIR + "/Python.h").st_mtime \
    or time_so < os.stat(os.path.abspath(__file__)).st_mtime:
    return True
  if os.path.isfile(main_path() + "/build/cmake_cython.cfg") \
    and time_so < os.stat(main_path() + "/build/cmake_cython.cfg").st_mtime:
    return True
  for file_end in ["pxd", "pyx", "cfg"]:
    time_in = os.stat(main_path() + "/cython/" + cython_from_cpp(constructor.global_loop) + "." \
      + file_end).st_mtime
    if time_so < time_in:
      return True
  dependent_files = [ x for x in constructor.include_files if "HyperHDG/" in x ]
  return check_dependent_files(dependent_files, time_so)

## \brief   Check whether any dependent files have been changed.
def check_dependent_files(dependent_files, time_so):
  for x in dependent_files:
    filename = main_path() + "/include/" + x
    if not os.path.isfile(filename):
      filename = main_path + x
    if time_so < os.stat(filename).st_mtime:
      return True
    with open(filename, 'r') as file:
      content = file.read()
    matches = re.finditer("\#include \<HyperHDG\/", content)
    for match in matches:
      start = match.start() + 10
      end   = start + re.search('>', content[start:len(content)]).end() - 1
      name  = content[start:end]
      if name not in dependent_files:
        dependent_files.append(name)
  return False

## Return path to main directory of HyperHDG.
def main_path():
  return os.path.dirname(os.path.abspath(__file__))

## Transform C++ class name to cython file name.
def cython_from_cpp(name):
  helper = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
  return re.sub('([a-z0-9])([A-Z])', r'\1_\2', helper).lower()