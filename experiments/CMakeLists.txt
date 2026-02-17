set(CMAKE_CXX_STANDARD 23)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter -Wno-format-security -Wno-format-truncation)

#######################

link_libraries(PkgConfig::PETSC)
find_package(HDF5)
if (HDF5_FOUND)
  link_libraries(HDF5::HDF5)
endif()

add_library(prin2 prin2.cxx)
add_library(net2as net2as.cxx petsc_parhip.cxx)
target_link_libraries(net2as PRIVATE prin2 parhip_interface_static)

add_executable(graph_gen graph_gen.cxx)
add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo nanoflann::nanoflann)

link_libraries(HyperHDG prin2)

add_executable(heat heat.cxx)
add_executable(wave wave.cxx)
add_executable(network network.cxx)

target_link_libraries(network net2as)
