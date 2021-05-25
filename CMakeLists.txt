cmake_minimum_required(VERSION 3.16)

project(HyperHDG LANGUAGES CXX)

# Include the submodule projects
add_subdirectory(submodules)

## Check whether appropriate compiler version is used
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient gcc/g++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient clang/clang++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
endif()

## Find packages that might become necessary.
find_package (LAPACK)

# Create a target that transports all the necessary information to downstream projects
add_library(HyperHDG INTERFACE)

# Add the correct include directories
target_include_directories(
  HyperHDG
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Link our target against upstream dependencies
target_link_libraries(
  HyperHDG
  INTERFACE
  HyperHDG::tpp
  HyperHDG::tpcc
  lapack
)

# Add the necessary C++ standard requirement
target_compile_features(
  HyperHDG
  INTERFACE
  cxx_std_20
)

# Add the examples subdirectory
add_subdirectory(examples)

# Add tests if requested through the BUILD_TESTING option (default: ON)
include(CTest)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests_c++)
endif()

# Build documentation if requested through the BUILD_DOC option (default: ON)
if(BUILD_DOC)
  add_subdirectory(doxygen)
endif()

# Add a namespaced alias target for the export
add_library(HyperHDG::HyperHDG ALIAS HyperHDG)

# Install targets, headers and configuration
include(GNUInstallDirs)
install(
  TARGETS HyperHDG
  EXPORT HyperHDG-config
)

install(
  EXPORT HyperHDG-config
  NAMESPACE HyperHDG::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HyperHDG
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

find_package (Python3 COMPONENTS Interpreter Development NumPy)
find_program (CYTHON_EXECUTABLE NAMES cython cython.bat cython3)

## Prepare nice output for include string.
get_directory_property (inc_dirs INCLUDE_DIRECTORIES)
list(TRANSFORM inc_dirs PREPEND "-I")
list(JOIN inc_dirs " " INCLUDE_STRING)

## Prepare nice output for linked libraries.
foreach (linked_lib ${LAPACK_LIBRARIES})
  string(APPEND LINKED_STRING ${linked_lib} " ")
endforeach()

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
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "libs    = ${LINKED_STRING}\n")

## Write commands for cython.

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "\n[cython]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "command = ${CYTHON_EXECUTABLE}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "flags   = -3 --cplus\n")

## Write information about python.

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "\n[python]\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "ver_maj = ${Python3_VERSION_MAJOR}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "ver_min = ${Python3_VERSION_MINOR}\n")
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg "dir     = ${Python3_INCLUDE_DIRS}\n")

if(NOT NOPYTHONTESTS)
  add_subdirectory(tests_python)
endif()

# Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
