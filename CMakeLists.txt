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
# On macOS prefer Apple's Accelerate framework, which bundles BLAS/LAPACK, so no external
# (e.g. keg-only Homebrew) libraries are required. On Linux this branch is skipped and the
# usual system LAPACK is discovered as before.
if(APPLE AND NOT DEFINED BLA_VENDOR)
  set(BLA_VENDOR Apple)
endif()
find_package (LAPACK REQUIRED)

## HDF5 for .geo.h5 domain files (include/HyperHDG/read_domain_hdf5.hxx). CMake's FindHDF5 looks
## at HDF5_ROOT and CMAKE_PREFIX_PATH (a conda environment: $CONDA_PREFIX; a spack environment sets
## it on activation) and at h5cc / h5pcc on PATH. A parallel HDF5 (e.g. hdf5-openmpi) includes
## mpi.h, so MPI is looked up as well in that case. The Cython modules (cmake_cython.cfg.in) are
## built without HDF5 and refuse .geo.h5 files with a message; use the nanobind path for them.
option(HYPERHDG_HDF5 "Read and write .geo.h5 domain files with HDF5 (if found)" ON)
if(HYPERHDG_HDF5)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    set(HYPERHDG_HDF5_FOUND TRUE)
    set(HYPERHDG_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    set(HYPERHDG_HDF5_LIBRARIES ${HDF5_C_LIBRARIES})
    if(HDF5_IS_PARALLEL)
      find_package(MPI REQUIRED COMPONENTS C)
      list(APPEND HYPERHDG_HDF5_INCLUDE_DIRS ${MPI_C_INCLUDE_DIRS})
      list(APPEND HYPERHDG_HDF5_LIBRARIES ${MPI_C_LIBRARIES})
    endif()
  else()
    message(WARNING "HDF5 not found: .geo.h5 domain files cannot be read. Set HDF5_ROOT (or "
                    "CMAKE_PREFIX_PATH) to its install prefix, or -DHYPERHDG_HDF5=OFF.")
  endif()
endif()
# On Windows, FindPython looks up interpreters registered in the Windows registry before the ones
# on PATH. This selects e.g. a python.org Python without NumPy instead of the MSYS2 Python of the
# shell, so prefer the interpreter on PATH, which is the one that runs the Python scripts.
if(WIN32 AND NOT DEFINED Python3_FIND_REGISTRY)
  set(Python3_FIND_REGISTRY LAST)
endif()
find_package (Python3 COMPONENTS Interpreter Development NumPy)

# Locate Cython. Searching $PATH alone misses installations inside a virtual environment or a
# Homebrew Python, so additionally hint the directory of the interpreter CMake selected, and fall
# back to invoking it as a module ("python -m cython") if no executable is found.
get_filename_component (_python_bindir "${Python3_EXECUTABLE}" DIRECTORY)
find_program (CYTHON_EXECUTABLE NAMES cython cython3 cython.bat HINTS "${_python_bindir}")
# On Windows, found "cython" scripts cannot necessarily be executed by the Windows shell. Invoking
# Cython as a module of the selected interpreter always works.
if(NOT CYTHON_EXECUTABLE OR WIN32)
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
  ${HYPERHDG_HDF5_INCLUDE_DIRS}
)
if(HYPERHDG_HDF5_FOUND)
  target_compile_definitions(HyperHDG INTERFACE HYPERHDG_HAVE_HDF5)
  target_compile_options(HyperHDG INTERFACE ${HDF5_C_DEFINITIONS})
endif()

## Link our target against upstream dependencies
target_link_libraries(
  HyperHDG
  INTERFACE
  HyperHDG::tpp
  HyperHDG::tpcc
  ${LAPACK_LIBRARIES}
  ${HYPERHDG_HDF5_LIBRARIES}
)

## Add the necessary C++ standard requirement
target_compile_features(
  HyperHDG
  INTERFACE
  cxx_std_20
)

## Optionally parallelize the assembly of matrices with OpenMP. On macOS,
## Apple Clang needs the libomp library, e.g. -DOpenMP_ROOT=$(brew --prefix libomp).
option(HYPERHDG_OPENMP "Parallelize matrix assembly with OpenMP" OFF)
if(HYPERHDG_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(HyperHDG INTERFACE OpenMP::OpenMP_CXX)
  # The Python interface compiles its modules with CMAKE_CXX_FLAGS (cf. cmake_cython.cfg.in).
  string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
  foreach(omp_include_dir ${OpenMP_CXX_INCLUDE_DIRS})
    string(APPEND CMAKE_CXX_FLAGS " -I${omp_include_dir}")
  endforeach()
  set(HYPERHDG_OPENMP_LIBRARIES ${OpenMP_CXX_LIBRARIES})
  # Let the Python modules find the OpenMP runtime (e.g. libomp outside the system directories).
  if(NOT WIN32)
    foreach(omp_lib ${OpenMP_CXX_LIBRARIES})
      get_filename_component(omp_lib_dir "${omp_lib}" DIRECTORY)
      string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,-rpath,${omp_lib_dir}")
    endforeach()
  endif()
endif()

## Add the examples subdirectory
add_subdirectory(examples)

## Add tests if requested through the BUILD_TESTING option (default: ON)
include(CTest)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests_c++)
endif()

## Add the command-line tools (geo2h5; its tests need enable_testing above)
add_subdirectory(tools)

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

## Python extension modules must not be linked against libpython on Linux and macOS: If the
## interpreter itself contains a static copy of libpython (e.g. Python managed by uv/pyenv, conda,
## python-build-standalone), a second Python runtime would be loaded and the import fails with
## "initialization ... did not return an extension module". The symbols are provided by the running
## interpreter instead (macOS requires an explicit permission for undefined symbols). On Windows,
## however, extension modules have to be linked against the Python import library.
if(WIN32)
  set(PYTHON_MODULE_LINK ${Python3_LIBRARIES})
elseif(APPLE)
  set(PYTHON_MODULE_LINK "-undefined dynamic_lookup")
endif()

## Prepare nice output for linked libraries.
foreach (linked_lib ${LAPACK_LIBRARIES} ${PYTHON_MODULE_LINK} ${HYPERHDG_OPENMP_LIBRARIES})
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
