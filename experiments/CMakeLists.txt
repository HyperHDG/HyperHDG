set(WARNINGS -Wall -Wextra -Werror -Wshadow -Wformat=2 -Wnull-dereference -Wuninitialized -Wmissing-declarations -Wno-unused-parameter -Wno-shadow -Wno-unused-but-set-variable  -Wno-unused-variable)

add_library(geobin geobin.cxx)
target_link_libraries(geobin PRIVATE bxzstr PRIVATE fmtlog-static)
target_compile_options(geobin PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(make_geo make_geo.cxx)
target_link_libraries(make_geo PRIVATE geobin PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static)
target_compile_options(make_geo PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>
)

add_executable(partitioner partitioner.cxx)
target_link_libraries(partitioner geobin KaHIP libkaffpa)
target_compile_options(partitioner PRIVATE ${WARNINGS} -std=gnu++2b $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)

add_executable(read_domain read_domain.cxx)
target_link_libraries(read_domain PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static)
target_compile_options(read_domain PRIVATE -std=c++23 ${WARNINGS} $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
target_include_directories(read_domain PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
