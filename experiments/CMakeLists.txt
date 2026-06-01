set(CMAKE_CXX_STANDARD 23)

find_package(HDF5)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter -Wno-format-security -Wno-format-truncation)

option(HYPERHDG_PARHIP "build HyperHDG with ParHIP" OFF)
set(PARHIP_ROOT "" CACHE PATH "ParHIP install prefix")

#######################

add_library(prin2 prin2.cxx)
target_link_libraries(prin2 PUBLIC PkgConfig::PETSC)

add_library(net2as net2as.cxx)
target_link_libraries(net2as
  PRIVATE
  prin2
  $<$<BOOL:${HDF5_FOUND}>:HDF5::HDF5>
  PUBLIC
  PkgConfig::PETSC
)

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

link_libraries(HyperHDG prin2)

add_executable(heat heat.cxx)
add_executable(wave wave.cxx)

link_libraries(net2as)

add_executable(timowave timowave.cxx)
add_executable(network network.cxx)
