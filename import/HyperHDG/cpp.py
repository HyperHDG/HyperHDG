import re

## \brief   Extract classname from name that might contain template arguemnts.
def extract_classname(fullname):
  index = re.search('\<|\ ', fullname)
  return fullname[0:index.start() if index != None else len(fullname)]

## \brief   Find file with definition of classname in directory folder.
def find_definition(folder, classname):
  for file in glob.glob(main_path() + "/include/HyperHDG/" + folder + "/*.hxx"):
    with open(file, "r") as hxxfile:
      if ("\nclass " + classname + "\n{" or "\nstruct " + classname + "\n{") in hxxfile.read():
        return re.sub(main_path() + "/include/", '', file)
  assert False, "File containing defintion of " + classname + " has not been found!"
