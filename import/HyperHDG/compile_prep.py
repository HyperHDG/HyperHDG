import glob, os, re
from .cmake import options
from .names import cython_from_cpp
from .paths import main_dir, this_dir

## \brief   Generate shell commands from CMAKE parameters.
def compile_commands(python_class, opt):
  assert isinstance(opt, options)
  cython_command = "cd " + main_dir() + "/build/cython_files/; " \
    + opt.cython_com + " " + opt.cython_flg + " " + python_class + ".pyx";
  compile_command = "cd " + main_dir() + "; " + opt.compile_com + " " + opt.compile_inc + " " + \
    opt.compile_flg + " --std=gnu++" + opt.compile_std + " -c ./build/cython_files/" + \
    python_class + ".cpp -o ./build/cython_files/" + python_class + ".o";
  link_command = "cd " + main_dir() + "; " \
    + opt.link_com + " " + opt.link_flg + " ./build/cython_files/" + python_class + ".o " \
    + " -o build/shared_objects/" + python_class + ".so " + opt.link_lib;
  return cython_command, compile_command, link_command

## \brief   Check whether recompilation of executable is necessary.
def need_compile(conf, python_class, opt):
  assert isinstance(opt, options)
  if own_code(conf, python_class):
    return True
  if not os.path.isfile(main_dir() + "/build/shared_objects/" + python_class + ".so"):
    return True
  time_so = os.stat(main_dir() + "/build/shared_objects/" + python_class + ".so").st_mtime
  assert os.path.isfile(opt.py_dir + "/Python.h"), "Include file for Python needs to exist!"
  if time_so < os.stat(opt.py_dir + "/Python.h").st_mtime \
    or time_so < os.stat(os.path.abspath(__file__)).st_mtime:
    return True
  if os.path.isfile(main_dir() + "/build/cmake_cython.cfg") \
    and time_so < os.stat(main_dir() + "/build/cmake_cython.cfg").st_mtime:
    return True
  for file_end in ["pxd", "pyx", "cfg"]:
    time_in = os.stat(main_dir() + "/cython/" + cython_from_cpp(conf.global_loop) + "." \
      + file_end).st_mtime
    if time_so < time_in:
      return True
  for file in glob.glob(this_dir()):
    if time_so < os.stat(file).st_mtime:
      return True
  dependent_files = [ x for x in conf.include_files if "HyperHDG/" in x ]
  return need_compile_check_hy_files(dependent_files, time_so)

## \brief   Check whether any dependent files have been changed.
def need_compile_check_hy_files(dependent_files, time_so):
  for x in dependent_files:
    filename = main_dir() + "/include/" + x
    if not os.path.isfile(filename):
      filename = main_dir + x
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

## \brief   Check whether own code differs from last iteration.
def own_code(conf, python_class):
  if conf.cpp_code == "":
    return False
  if os.path.isfile(main_dir() + "/build/cython_files/" + python_class + ".hxx"):
    with  open(main_dir() + "/build/cython_files/" + python_class + ".hxx", "r") as file:
      content = file.read()
    if content == conf.cpp_code:
      return False
  with  open(main_dir() + "/build/cython_files/" + python_class + ".hxx", "w") as file:
    file.write(conf.cpp_code)
  return True
