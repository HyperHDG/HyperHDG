import configparser, os, sys
from .paths import main_dir

## \brief   Object that comprises all information for HyperHDG to create a problem.
class options:
  compile_com = ""
  compile_inc = ""
  compile_flg = ""
  compile_std = ""
  link_com    = ""
  link_flg    = ""
  link_lib    = ""
  cython_com  = ""
  cython_flg  = ""
  py_ver_maj  = ""
  py_ver_min  = ""
  py_dir      = ""

## \brief   Read out parameters for Cython compilation provided by CMAKE.
def get_options():
  opt = options()
  if not os.path.isfile(main_dir() + "/build/cmake_cython.cfg") \
   and not os.path.isfile("/home/runner/work/HyperHDG/build/cmake_cython.cfg"):
    print("CMAKE files do not exist ... using default values for Cython!")
    opt.compile_com = "g++-10"
    opt.compile_inc = "-I. -Iinclude -Isubmodules/tensor_product_chain_complex.git/include \
      -Isubmodules/tensor_product_polynomials.git/include \
      -I/usr/include/python" + str(sys.version_info.major) + "." + str(sys.version_info.minor)
    opt.compile_flg = "-pthread -g -fwrapv -O2 -Wall -pedantic -g -fstack-protector-strong \
      -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC"
    opt.compile_std = "20"
    opt.link_com = "g++-10"
    opt.link_flg = "-pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
      -Wl,-z,relro -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong -Wformat \
      -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2"
    opt.link_lib = "-llapack"
    opt.cython_com = "cython3"
    opt.cython_flg = "-3 --cplus"
    opt.py_ver_maj = str(sys.version_info.major)
    opt.py_ver_min = str(sys.version_info.minor)
    opt.py_dir = "/usr/include/python"+str(sys.version_info.major)+"."+str(sys.version_info.minor)
  else:
    config = configparser.ConfigParser()
    if os.path.isfile(main_dir() + "/build/cmake_cython.cfg"):
      config.read(main_dir() + "/build/cmake_cython.cfg")
    else:  # os.path.isfile("/home/runner/work/HyperHDG/build/cmake_cython.cfg")
      config.read("/home/runner/work/HyperHDG/build/cmake_cython.cfg")
    opt.compile_com = config['compiler']['command']
    opt.compile_inc = config['compiler']['includes']
    opt.compile_flg = config['compiler']['flags']
    opt.compile_std = config['compiler']['standard']
    opt.link_com = config['linker']['command']
    opt.link_flg = config['linker']['flags']
    opt.link_lib = config['linker']['libs']
    opt.cython_com = config['cython']['command']
    opt.cython_flg = config['cython']['flags']
    opt.py_ver_maj = config['python']['ver_maj']
    opt.py_ver_min = config['python']['ver_min']
    opt.py_dir = config['python']['dir']
  assert int(opt.py_ver_maj), "Python versions below 3 are not supported!"
  assert opt.py_dir in opt.compile_inc, "Python directory must be included!"
  assert int(opt.py_ver_maj) == sys.version_info.major and \
    int(opt.py_ver_min) == sys.version_info.minor, "Python version is not CMAKE's Python version!"
  return opt
