set(WARNINGS -Wall -Wextra -Werror -Wshadow -Wformat=2 -Wnull-dereference -Wuninitialized -Wmissing-declarations -Wfloat-equal -Wno-unused-parameter -Wno-shadow)

add_executable(make_geo make_geo.cxx)
target_compile_features(make_geo PUBLIC cxx_std_23)
target_link_libraries(make_geo PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static PRIVATE bxzstr)
target_compile_options(make_geo PRIVATE ${WARNINGS} $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>
)

add_executable(read_domain read_domain.cxx)
target_link_libraries(read_domain PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static)
target_compile_options(read_domain PRIVATE -std=c++23 ${WARNINGS} $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>)
target_include_directories(read_domain PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
