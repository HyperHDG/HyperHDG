import os, sys, importlib, re, datetime
from .cmake import get_options
from .compile_prep import compile_commands, need_compile
from .config import *
from .names import files, cython_from_cpp
from .paths import main_dir

## \brief   Function to import classes of the HyperHDG package using Cython.
#
#  \param   conf hyperhdg_conf object.
#
#  \authors   Guido Kanschat, Heidelberg University, 2020.
#  \authors   Andreas Rupp, Heidelberg University, 2020.
def include(conf):
  start_time = datetime.datetime.now()
  # Check that conf is appropriatelly filled.
  assert isinstance(conf, config) and consistent(conf)

  # Start program.
  options = get_options()
  print("Cythonizing ... ", end='', flush=True)

  # Create folders and log files and check for consistency.
  os.system("mkdir -p " + main_dir() + "/build " + main_dir() + "/build/cython_files " \
            + main_dir() + "/build/shared_objects")
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
      content = re.sub("C\+\+ClassName", "\"" + cpp_class + "\"", content)
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
      compile_command += " -DNDEBUG";
    if not (conf.allow_file_output):
      compile_command += " -DNOFILEOUT"
    #Actually compile the prepared files.
    assert os.system(cython_command) == 0
    assert os.system(compile_command) == 0
    assert os.system(link_command) == 0

  try:
    mod = importlib.import_module(python_class)
  except ImportError as error:
    sys.path.append(main_dir() + "/build/shared_objects")
    mod = importlib.import_module(python_class)

  delta_time_ms = 1000 * (datetime.datetime.now() - start_time).total_seconds()
  if compilation_necessary:
    print("DONE with compilation in " + "{:.2f}".format(delta_time_ms) + " milliseconds.")
  else:
    print("DONE without compilation in " + "{:.2f}".format(delta_time_ms) + " milliseconds.")

  return getattr(mod, python_class)
