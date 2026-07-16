cmake_minimum_required(VERSION 3.16)

option(HYPERHDG_PETSC "build with PETSc" OFF)
option(HYPERHDG_MPI "build with MPI" OFF)
option(HYPERHDG_PARHIP "build HyperHDG with ParHIP" OFF)
set(PARHIP_ROOT "" CACHE PATH "ParHIP install prefix")

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
  add_subdirectory(experiments)
  find_package(HDF5)
else()
  message("NO PETSC - no experiments")
endif()

if(HYPERHDG_MPI)
  find_package(MPI REQUIRED)
else()
  message("NO MPI")
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
find_package (LAPACK REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development NumPy)
find_program (CYTHON_EXECUTABLE NAMES cython cython.bat cython3)

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
  $<$<BOOL:${HYPERHDG_MPI}>:MPI::MPI_CXX>
  $<$<BOOL:${HDF5_FOUND}>:HDF5::HDF5>
)

if (HYPERHDG_PARHIP)
  # Sequential KaHIP (kaffpa) used by distribute_domain.hxx to partition the network graph.
  # Same install prefix as ParHIP (spack provides both). Linked into the header-only HyperHDG
  # target so every consumer that includes read_domain.hxx gets the symbol + header + define.
  find_library(KAHIP_LIBRARY NAMES kahip_static HINTS ${PARHIP_ROOT}/lib REQUIRED)
  find_path(KAHIP_INCLUDE_DIR kaHIP_interface.h HINTS ${PARHIP_ROOT}/include REQUIRED)
  message("FOUND KAHIP_LIBRARY=${KAHIP_LIBRARY}")
  # KaHIP is a build-time dependency only: HyperHDG is a header-only, installed/exported INTERFACE
  # target, so its public install interface must not carry the (here source-tree-prefixed) kahip
  # paths. $<BUILD_INTERFACE:...> applies these to in-tree consumers (our executables) but drops
  # them from the install export.
  target_link_libraries(HyperHDG INTERFACE $<BUILD_INTERFACE:${KAHIP_LIBRARY}>)
  target_include_directories(HyperHDG INTERFACE $<BUILD_INTERFACE:${KAHIP_INCLUDE_DIR}>)
  target_compile_definitions(HyperHDG INTERFACE $<BUILD_INTERFACE:HYPERHDG_PARHIP>)
endif()

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

## Prepare nice output for linked libraries.
foreach (linked_lib ${LAPACK_LIBRARIES} ${Python3_LIBRARIES})
  string(APPEND LINKED_STRING ${linked_lib} " ")
endforeach()

## Write a configuration file for Cython
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake_cython.cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg
)

if(NOT NOPYTHONTESTS)
  add_subdirectory(tests_python)
endif()

## Explicit nanobind python modules (petsc-free, see python/CMakeLists.txt)
add_subdirectory(python)

## Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
