set(WARNINGS -Wall -Wextra -Werror -Wshadow -Wformat=2 -Wnull-dereference -Wuninitialized -Wmissing-declarations -Wno-unused-parameter -Wno-shadow -Wno-unused-but-set-variable  -Wno-unused-variable)

add_library(stats stats.cxx)
target_compile_options(stats PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_library(geobin geobin.cxx)
target_link_libraries(geobin PRIVATE bxzstr PRIVATE fmtlog PUBLIC KaHIP)
target_compile_options(geobin PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo PRIVATE geobin PRIVATE nanoflann::nanoflann PRIVATE fmtlog)
target_compile_options(make_geo PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>
)

add_executable(partitioner partitioner.cxx)
target_link_libraries(partitioner geobin stats KaHIP libkaffpa CLI11::CLI11 bxzstr METIS::METIS fmtlog)
target_compile_options(partitioner PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(graph_gen graph_gen.cxx)
target_link_libraries(graph_gen CLI11::CLI11 KaHIP geobin fmtlog)
target_compile_options(graph_gen PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(read_domain read_domain.cxx)
target_link_libraries(read_domain nanoflann::nanoflann fmtlog geobin)
target_compile_options(read_domain PRIVATE -std=c++23 ${WARNINGS} $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
target_include_directories(read_domain PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
