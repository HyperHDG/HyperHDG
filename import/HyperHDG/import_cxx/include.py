import os, sys, importlib, re, datetime
from .cmake import get_options
from .compile_prep import compile_commands, need_compile, run_command
from .config import *
from .names import files, cython_from_cpp
from .paths import main_dir

## \brief   Function to import classes of the HyperHDG package using Cython.
#
#  \param   conf  Object of type config.
def include(conf):
  start_time = datetime.datetime.now()
  # Check that conf is appropriatelly filled.
  assert isinstance(conf, config) and consistent(conf)

  # Start program.
  options = get_options()
  print("Cythonizing ... ", end='', flush=True)

  # Create folders and log files and check for consistency.
  os.makedirs(main_dir() + "/build/cython_files", exist_ok=True)
  os.makedirs(main_dir() + "/build/shared_objects", exist_ok=True)
  if not os.path.isfile(main_dir() + "/build/cython_log.txt"):
    file = open(main_dir() + "/build/cython_log.txt", "w")
    file.write("Python version: " + str(sys.version_info))
    file.close()
  else:
    file = open(main_dir() + "/build/cython_log.txt", "r")
    assert file.readline() == "Python version: " + str(sys.version_info), \
           "Wrong Python version!\n \
            Remove build directory and reinstall HyperHDG or use previous Python version."
    file.close()

  # Evaluate configfile, define include files, and define dependent files.
  cy_replace = generate_cy_replace(conf)
  cpp_class, python_class, cython_class = files(conf)
  if conf.cpp_code == "":
    include_string = extract_includes(conf)
  else:
    include_string = extract_includes(conf, main_dir()+"/build/cython_files/"+python_class+".hxx")

  compilation_necessary = need_compile(conf, python_class, options)
  if compilation_necessary:
    # Copy pyx and pxd files from cython directory.
    for file_end in ["pxd", "pyx"]:
      with open(main_dir() + "/cython/" + cython_from_cpp(conf.global_loop) + "." \
        + file_end, "r") as file:
        content = file.read()
      content = re.sub(r"C\+\+ClassName", "\"" + cpp_class + "\"", content)
      content = re.sub("CythonClassName", cython_class, content)
      content = re.sub("PythonClassName", python_class, content)
      content = re.sub("IncludeFiles", include_string, content)
      for i in range(len(cy_replace)):
        content = re.sub("CyReplace" + '%02d' % (i+1), cy_replace[i], content)
      with open(main_dir() + "/build/cython_files/" + python_class + "." + file_end, "w") as file:
        file.write(content)
    # Prepare the compilation commands.
    cython_command, compile_command, link_command = compile_commands(python_class, options)
    if not (conf.debug_mode):
      compile_command += " -DNDEBUG"
    #Actually compile the prepared files.
    assert run_command(cython_command) == 0
    assert run_command(compile_command) == 0
    assert run_command(link_command) == 0

  # Since Python 3.8, Windows does not search PATH for the DLLs an extension module depends on
  # (libstdc++, LAPACK, ...). They reside next to the compiler (e.g. MSYS2's ucrt64/bin).
  if os.name == "nt" and hasattr(os, "add_dll_directory"):
    compiler_dir = os.path.dirname(options.compile_com)
    if os.path.isdir(compiler_dir):
      os.add_dll_directory(compiler_dir)

  try:
    mod = importlib.import_module(python_class)
  except ImportError as error:
    sys.path.append(main_dir() + "/build/shared_objects")
    mod = importlib.import_module(python_class)

  delta_t = 1000 * (datetime.datetime.now() - start_time).total_seconds()
  if compilation_necessary:
    print("DONE with compilation in " + "{:,.2f}".format(delta_t) + " milliseconds.", flush=True)
  else:
    print("DONE without compilation in " + "{:,.2f}".format(delta_t) + " milliseconds.", flush=True)

  return getattr(mod, python_class)
