set(CMAKE_CXX_STANDARD 23)

# apply options privately to all following libs/exe
add_compile_options(-Wall -Wextra -Werror -Wno-unused-parameter)

########## LIBS ############

add_library(libpartition libpartition.cxx)
target_include_directories(libpartition PUBLIC .)
target_link_libraries(libpartition PUBLIC kahip PUBLIC METIS_METIS)

add_library(stats stats.cxx)

add_library(geobin geobin.cxx)

########## EXE #############

add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo geobin nanoflann::nanoflann)

add_executable(partitioner partitioner.cxx)
target_link_libraries(partitioner geobin stats CLI11::CLI11 libpartition spdlog::spdlog)
target_compile_options(partitioner PRIVATE -Wno-deprecated-literal-operator)

add_executable(graph_gen graph_gen.cxx)
target_link_libraries(graph_gen geobin)

link_libraries(HyperHDG PkgConfig::PETSC)
add_executable(heat heat.cxx)
add_executable(wave wave.cxx)
add_executable(wave_conv wave_conv.cxx)
add_executable(network network.cxx)
target_link_libraries(network geobin)
