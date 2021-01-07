import os

## Return path to main directory of HyperHDG.
def main_dir():
  return os.path.dirname(os.path.abspath(__file__)) + "/../.."
