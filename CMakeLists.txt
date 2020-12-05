cmake_minimum_required(VERSION 3.16)
enable_testing()
set(default_build_type "Debug")

project(HyperHDG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic -Wextra -D_GLIBCXX_DEBUG")

find_package (LAPACK)
find_package (Python3 COMPONENTS Interpreter Development NumPy)
find_program( CYTHON_EXECUTABLE
  NAMES cython cython.bat cython3)

include_directories(include submodules/tensor_product_chain_complex.git/include/)

# This 
get_directory_property (inc_dirs INCLUDE_DIRECTORIES)
list(TRANSFORM inc_dirs PREPEND "-I ")
list(JOIN inc_dirs " " INCLUDE_STRING)

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

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output ${CMAKE_CURRENT_BINARY_DIR}/cython_files ${CMAKE_CURRENT_BINARY_DIR}/shared_objects)
file(COPY domains reproducables_python DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(tests_c++)
# add_subdirectory(doc)