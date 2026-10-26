import os

## Return path to directory of file.
def this_dir():
  return os.path.dirname(os.path.abspath(__file__))

## Return path to main directory of HyperHDG.
def main_dir():
  return os.path.abspath(this_dir() + "/../../..").replace(os.sep, "/")

## Return file ending of Python extension modules ('.pyd' on Windows, '.so' otherwise).
def so_suffix():
  return ".pyd" if os.name == "nt" else ".so"
