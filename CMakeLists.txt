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
find_package (LAPACK REQUIRED)
find_package (OpenACC)

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

find_package (Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_program (CYTHON_EXECUTABLE NAMES cython cython.bat cython3 REQUIRED)

## Add OpenMP compiler and linker flags.
foreach (linker_flag ${OpenACC_CXX_FLAGS})
  string(APPEND CMAKE_SHARED_LINKER_FLAGS ${linker_flag} " ")
  string(APPEND CMAKE_CXX_FLAGS ${linker_flag} " ")
endforeach()

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

# Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
