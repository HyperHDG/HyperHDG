set(CMAKE_CXX_STANDARD 23)

find_package(HDF5)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter -Wno-format-security -Wno-format-truncation)

#######################

add_library(petsc_util petsc_util.cxx)
target_link_libraries(petsc_util PUBLIC PkgConfig::PETSC)

if(NOT HYPERHDG_COMPLEX)
# cholmod comes with petsc's suite-sparse, but petsc.pc only carries it in Libs.private
find_library(CHOLMOD_LIBRARY cholmod HINTS ${PETSC_STATIC_LIBRARY_DIRS} REQUIRED)
cmake_path(GET CHOLMOD_LIBRARY PARENT_PATH CHOLMOD_LIBDIR)
find_path(CHOLMOD_INCLUDE_DIR cholmod.h HINTS ${CHOLMOD_LIBDIR}/../include REQUIRED)
add_executable(frontal_distribution frontal_distribution.cxx)
target_link_libraries(frontal_distribution PkgConfig::PETSC ${CHOLMOD_LIBRARY})
target_include_directories(frontal_distribution PRIVATE ${CHOLMOD_INCLUDE_DIR})
add_library(net2as net2as.cxx)
target_link_libraries(net2as
  PRIVATE
  petsc_util
  $<$<BOOL:${HDF5_FOUND}>:HDF5::HDF5>
  PUBLIC
  PkgConfig::PETSC
)
endif()

if (HYPERHDG_PARHIP)
  message("USING PARHIP_ROOT=${PARHIP_ROOT}")
  find_library(PARHIP_LIBRARY parhip_interface HINTS ${PARHIP_ROOT}/lib REQUIRED)
  find_path(PARHIP_INCLUDE_DIR parhip_interface.h HINTS ${PARHIP_ROOT}/include REQUIRED)
  message("FOUND PARHIP_LIBRARY=${PARHIP_LIBRARY}")

  add_library(petsc_parhip petsc_parhip.cxx)
  target_link_libraries(petsc_parhip PRIVATE ${PARHIP_LIBRARY} PUBLIC PkgConfig::PETSC)
  target_include_directories(petsc_parhip PRIVATE ${PARHIP_INCLUDE_DIR})
  target_link_libraries(net2as PRIVATE petsc_parhip)
  target_compile_definitions(net2as PUBLIC HYPERHDG_PARHIP)
endif()

link_libraries(HyperHDG petsc_util)

# The complex-PETSc build (preset "complex") serves the multi-stage Gauss timowave only; the
# remaining binaries (and net2as) are real-scalar production code.
if(NOT HYPERHDG_COMPLEX)
  add_executable(heat heat.cxx)
  add_executable(wave wave.cxx)
  link_libraries(net2as)
endif()

add_executable(timowave timowave.cxx)

if(NOT HYPERHDG_COMPLEX)
  add_executable(network network.cxx)
endif()
