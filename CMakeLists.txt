cmake_minimum_required(VERSION 3.16)
enable_testing()
set(default_build_type "Debug")

project(HyperHDG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Set the flags for the compiler and the linker.

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
  -pthread -g -fwrapv -O2 -Wall -pedantic -Wextra -g \
  -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} \
  -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro \
  -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong -Wformat \
  -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic -Wextra -D_GLIBCXX_DEBUG")

## Find packages that might become necessary.

find_package (LAPACK)
find_package (Python3 COMPONENTS Interpreter Development NumPy)
find_program (CYTHON_EXECUTABLE NAMES cython cython.bat cython3)

## Set directories to be included.

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/submodules/tensor_product_chain_complex.git/include/
  ${Python3_INCLUDE_DIRS})

# CONFIGURE_FILE(
#   ${CMAKE_CURRENT_SOURCE_DIR}/cython_import.py.in
#   ${CMAKE_CURRENT_BINARY_DIR}/cython_import.py
#   @ONLY)

# CONFIGURE_FILE(
#   ${CMAKE_CURRENT_SOURCE_DIR}/cython/cythonize.cxx.in
#   ${CMAKE_CURRENT_BINARY_DIR}/cython/cythonize.cxx
#   @ONLY)

# add_executable(cythonize ${CMAKE_CURRENT_BINARY_DIR}/cython/cythonize.cxx)
# target_compile_definitions(cythonize PUBLIC PYVERMAJ=${Python3_VERSION_MAJOR} PYVERMIN=${Python3_VERSION_MINOR})

# file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output ${CMAKE_CURRENT_BINARY_DIR}/cython_files ${CMAKE_CURRENT_BINARY_DIR}/shared_objects)
# file(COPY domains reproducables_python DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


## Prepare nice output for include string.

get_directory_property (inc_dirs INCLUDE_DIRECTORIES)
list(TRANSFORM inc_dirs PREPEND "-I")
list(JOIN inc_dirs " " INCLUDE_STRING)

## Write commands for compiler.

file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "[compiler]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "command  = ${CMAKE_CXX_COMPILER}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "includes = ${INCLUDE_STRING}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "flags    = ${CMAKE_CXX_FLAGS}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "standard = ${CMAKE_CXX_STANDARD}\n")

## Write commands for linker.

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "\n[linker]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "command = ${CMAKE_CXX_COMPILER}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "flags   = ${CMAKE_SHARED_LINKER_FLAGS}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "libs    = ${LAPACK_lapack_LIBRARY}\n")

## Write commands for cython.

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "\n[cython]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "command = ${CYTHON_EXECUTABLE}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "flags   = -3 --cplus\n")

## Write information about python.

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "\n[python]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "ver_maj = ${Python3_VERSION_MAJOR}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "ver_min = ${Python3_VERSION_MINOR}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "dir     = ${Python3_INCLUDE_DIRS}\n")


## Add subdirectories of tests.

add_subdirectory(tests_c++)
if(NOT DEFINED NOPYTHONTESTS)
  add_subdirectory(tests_python)
endif()
# add_subdirectory(doc)
