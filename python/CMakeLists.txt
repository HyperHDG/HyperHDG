# Explicit nanobind extension modules: a tiny .cxx instantiates bind_python<> (see
# include/HyperHDG/bind_python.hxx) for the template combinations actually needed and CMake
# builds a plain .so that python imports -- no runtime compile, no PETSc/MPI linkage.

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/submodules/nanobind.git/CMakeLists.txt)
  message(STATUS "nanobind submodule missing (run 'git submodule update --init --recursive') - skipping python modules")
  return()
endif()

find_package(Python 3.9 COMPONENTS Interpreter Development.Module)
if(NOT Python_FOUND)
  message(STATUS "no python Development.Module - skipping python modules")
  return()
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/nanobind.git nanobind)

nanobind_add_module(timowave_py timowave_py.cxx)
target_compile_features(timowave_py PRIVATE cxx_std_20)
# HyperHDG headers without the HyperHDG interface target: that target would drag in the
# PETSc/MPI link interface and the HYPERHDG_PETSC/HYPERHDG_MPI definitions -- this module is
# deliberately the serial, petsc-free compile of the very same headers.
target_include_directories(
  timowave_py
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/experiments
)
target_link_libraries(timowave_py PRIVATE HyperHDG::tpp HyperHDG::tpcc ${LAPACK_LIBRARIES})
