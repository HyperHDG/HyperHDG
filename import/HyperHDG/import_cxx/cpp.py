import glob, os, re
from .paths import main_dir

## \brief   Extract classname from name that might contain template arguemnts.
def extract_classname(fullname):
  index = re.search(r'\<|\ ', fullname)
  return fullname[0:index.start() if index != None else len(fullname)]

## \brief   Find file with definition of classname in directory folder.
def find_definition(folder, classname):
  # A definition (not a declaration) of the class, possibly with base classes.
  definition = re.compile(r"\n(class|struct) " + re.escape(classname) + r"\s*(:[^;{]*)?\{")
  for file in glob.glob(main_dir() + "/include/HyperHDG/" + folder + "/*.hxx"):
    with open(file, "r") as hxxfile:
      if definition.search(hxxfile.read()):
        # Path relative to the include directory using '/' (glob yields '\' on Windows).
        return os.path.relpath(file, main_dir() + "/include").replace(os.sep, "/")
  assert False, "File containing defintion of " + classname + " has not been found!"
