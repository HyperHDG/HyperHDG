# https://stackoverflow.com/questions/16993927/using-cython-to-link-python-to-a-shared-library

# setup.py file
import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# build "myext.so" python extension to be added to "PYTHONPATH" afterwards...
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize("ClassWrapper.pyx", language="c++")
#    ext_modules = [
#        Extension("ClassWrapper", 
#                  sources=["ClassWrapper.pyx"],
#                  language="c++",                   # remove this if C and not C++
#                  extra_compile_args=["--std=c++17"]
#             )
#        ]
)
