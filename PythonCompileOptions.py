# https://stackoverflow.com/questions/16993927/using-cython-to-link-python-to-a-shared-library

# setup.py file
import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# clean previous build
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if (name.startswith("ClassWrapper") and not(name.endswith(".pyx") or name.endswith(".pxd"))):
            os.remove(os.path.join(root, name))
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)

# build "myext.so" python extension to be added to "PYTHONPATH" afterwards...
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("ClassWrapper", 
                  sources=["ClassWrapper.pyx",
                           "HyperEdge.C",
                           "Topology.C",
                           "DiffusionProblem.C",
                           "DiffusionSolver.C",
                           "FuncAndQuad.C",
                           "HyperGraph.C",
                           "HyperNodeFactory.C",
                           "HyperEdge_Geometry.C",
                           "Geometry.C",
                           "Plotter.C",
                           "Point.C"
                       ],
                  language="c++",                   # remove this if C and not C++
                  extra_compile_args=["--std=c++17"],
                  extra_link_args=["-llapack"]
             )
        ]
)
