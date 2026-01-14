set(CMAKE_CXX_STANDARD 23)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter -Wno-format-security -Wno-format-truncation)

########## EXE #############

add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo nanoflann::nanoflann PkgConfig::PETSC)

add_executable(graph_gen graph_gen.cxx)
target_link_libraries(graph_gen PkgConfig::PETSC)

link_libraries(HyperHDG PkgConfig::PETSC)

add_executable(heat heat.cxx)
add_executable(wave wave.cxx)
add_executable(network network.cxx)
