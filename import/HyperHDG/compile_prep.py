## \brief   Read out parameters for Cython compilation provided by CMAKE.
def cmakes_options():
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
    CYTHON_COM = "cython"
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

## \brief   Generate shell commands from CMAKE parameters.
def compile_commands(python_class):
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
def need_compile(conf, python_class):
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
    time_in = os.stat(main_path() + "/cython/" + cython_from_cpp(conf.global_loop) + "." \
      + file_end).st_mtime
    if time_so < time_in:
      return True
  dependent_files = [ x for x in conf.include_files if "HyperHDG/" in x ]
  return need_compile_check_hy_files(dependent_files, time_so)

## \brief   Check whether any dependent files have been changed.
def need_compile_check_hy_files(dependent_files, time_so):
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
