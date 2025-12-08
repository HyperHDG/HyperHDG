set(WARNINGS -Wall -Wextra -Werror -Wshadow -Wformat=2 -Wnull-dereference -Wuninitialized -Wmissing-declarations -Wno-unused-parameter -Wno-shadow -Wno-unused-but-set-variable  -Wno-unused-variable -Wno-deprecated-literal-operator)
# TODO: use add_compile_options?

########## LIBS ############

add_library(libpartition libpartition.cxx)
target_compile_options(libpartition PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
target_include_directories(libpartition PUBLIC .)
target_link_libraries(libpartition PUBLIC spdlog PUBLIC kahip PUBLIC METIS_METIS)

add_library(stats stats.cxx)
target_compile_options(stats PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_library(geobin geobin.cxx)
target_link_libraries(geobin PUBLIC bxzstr)
target_compile_options(geobin PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
target_compile_definitions(geobin PUBLIC HAVE_GEOBIN)

########## EXE #############

add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo PRIVATE geobin PRIVATE nanoflann::nanoflann PRIVATE fmtlog)
target_compile_options(make_geo PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>
)

add_executable(partitioner partitioner.cxx)
target_link_libraries(partitioner geobin stats CLI11::CLI11 libpartition spdlog::spdlog)
target_compile_options(partitioner PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(graph_gen graph_gen.cxx)
target_link_libraries(graph_gen CLI11::CLI11 geobin fmtlog)
target_compile_options(graph_gen PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(heat heat.cxx)
target_link_libraries(heat PRIVATE HyperHDG PRIVATE geobin PRIVATE PkgConfig::PETSC)
target_compile_options(heat PRIVATE -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(wave wave.cxx)
target_link_libraries(wave PRIVATE HyperHDG PRIVATE geobin PRIVATE PkgConfig::PETSC)
target_compile_options(wave PRIVATE -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(wave_conv wave_conv.cxx)
target_link_libraries(wave_conv PRIVATE HyperHDG PRIVATE geobin PRIVATE PkgConfig::PETSC)
target_compile_options(wave_conv PRIVATE -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(network network.cxx)
target_link_libraries(network PRIVATE HyperHDG PRIVATE geobin PRIVATE PkgConfig::PETSC)
target_compile_options(network PRIVATE -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
