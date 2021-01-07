import os

## Return path to directory of file.
def this_dir():
  return os.path.dirname(os.path.abspath(__file__))

## Return path to main directory of HyperHDG.
def main_dir():
  return this_dir() + "/../.."
