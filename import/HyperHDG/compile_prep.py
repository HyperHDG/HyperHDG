

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
