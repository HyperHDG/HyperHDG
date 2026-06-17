cmake_minimum_required(VERSION 3.16)

option(HYPERHDG_PETSC "build with PETSc" OFF)
option(HYPERHDG_MPI "build with MPI" OFF)

# the compilers must be set before project() is called
if (HYPERHDG_PETSC)
  set(HYPERHDG_PETSC_PREFIX "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}" CACHE STRING "PETSc prefix")
  message("USING PETSC_PREFIX=${HYPERHDG_PETSC_PREFIX}")
  set(ENV{PKG_CONFIG_PATH} ${HYPERHDG_PETSC_PREFIX}/lib/pkgconfig)

  find_package(PkgConfig REQUIRED)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
  pkg_get_variable(CMAKE_C_COMPILER PETSc ccompiler)
  pkg_get_variable(CMAKE_CXX_COMPILER PETSc cxxcompiler)
endif()

project(HyperHDG LANGUAGES CXX)

# Include the submodule projects
add_subdirectory(submodules)
include_directories(submodules)

if(HYPERHDG_PETSC)
  add_subdirectory(reproducibles_petsc)
  find_package(HDF5)
else()
  message("--- NO PETSC")
endif()

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
# On macOS prefer Apple's Accelerate framework, which bundles BLAS/LAPACK, so no external
# (e.g. keg-only Homebrew) libraries are required. On Linux this branch is skipped and the
# usual system LAPACK is discovered as before.
if(APPLE AND NOT DEFINED BLA_VENDOR)
  set(BLA_VENDOR Apple)
endif()
find_package (LAPACK REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development NumPy)

# Locate Cython. Searching $PATH alone misses installations inside a virtual environment or a
# Homebrew Python, so additionally hint the directory of the interpreter CMake selected, and fall
# back to invoking it as a module ("python -m cython") if no executable is found.
get_filename_component (_python_bindir "${Python3_EXECUTABLE}" DIRECTORY)
find_program (CYTHON_EXECUTABLE NAMES cython cython3 cython.bat HINTS "${_python_bindir}")
if(NOT CYTHON_EXECUTABLE)
  set(CYTHON_EXECUTABLE "${Python3_EXECUTABLE} -m cython")
endif()

## Create a target that transports all the necessary information to downstream projects
add_library(HyperHDG INTERFACE)

## Add the correct include directories
target_include_directories(
  HyperHDG
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

## Link our target against upstream dependencies
target_link_libraries(
  HyperHDG
  INTERFACE
  HyperHDG::tpp
  HyperHDG::tpcc
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${HYPERHDG_PETSC}>:PkgConfig::PETSC>
  $<$<BOOL:${HDF5_FOUND}>:HDF5::HDF5>
)

## Add the necessary C++ standard requirement
target_compile_features(
  HyperHDG
  INTERFACE
  cxx_std_20
)

target_compile_definitions(
  HyperHDG
  INTERFACE
  $<$<BOOL:${HYPERHDG_MPI}>:HYPERHDG_MPI>
  $<$<BOOL:${HYPERHDG_PETSC}>:HYPERHDG_PETSC>
)

## Add the examples subdirectory
add_subdirectory(examples)

## Add tests if requested through the BUILD_TESTING option (default: ON - set by CTest)
include(CTest)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests_c++)
endif()

## Build documentation if requested through the BUILD_DOC option (default: ON)
if(BUILD_DOC)
  add_subdirectory(doxygen)
endif()

## Add a namespaced alias target for the export
add_library(HyperHDG::HyperHDG ALIAS HyperHDG)

## Install targets, headers and configuration
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

## On macOS the Cython compile/link step invokes the compiler directly (see
## import/HyperHDG/import_cxx/compile_prep.py) rather than through CMake, so it does not inherit
## the SDK sysroot that CMake passes to its own targets. Without it, system headers such as
## <assert.h> are not found. Forward CMake's detected sysroot so the manual invocation matches.
## On Linux this variable stays empty and the generated config is unchanged.
if(APPLE AND CMAKE_OSX_SYSROOT)
  set(MACOS_SYSROOT_FLAG "-isysroot ${CMAKE_OSX_SYSROOT}")
endif()

## Prepare nice output for linked libraries.
foreach (linked_lib ${LAPACK_LIBRARIES} ${Python3_LIBRARIES})
  # On macOS, find_package(LAPACK) reports Accelerate as a path to the .framework directory, which
  # the linker cannot consume directly. Translate such entries into the proper "-framework <name>"
  # flag. Regular library paths (the Linux case) are passed through unchanged.
  if(linked_lib MATCHES "/([^/]+)\\.framework/?$")
    string(APPEND LINKED_STRING "-framework " ${CMAKE_MATCH_1} " ")
  else()
    string(APPEND LINKED_STRING ${linked_lib} " ")
  endif()
endforeach()

## Write a configuration file for Cython
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake_cython.cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg
)

if(NOT NOPYTHONTESTS)
  add_subdirectory(tests_python)
endif()

## Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
