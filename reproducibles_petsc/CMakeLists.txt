set(CMAKE_CXX_STANDARD 23)

find_package(HDF5)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter -Wno-format-security -Wno-format-truncation)

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

link_libraries(HyperHDG prin2)

add_executable(heat heat.cxx)
add_executable(wave wave.cxx)

link_libraries(net2as)

add_executable(timowave timowave.cxx)
add_executable(network network.cxx)
