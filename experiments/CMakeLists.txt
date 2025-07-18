set(WARNINGS -Wall -Wextra -Werror -Wshadow -Wformat=2 -Wnull-dereference -Wuninitialized -Wmissing-declarations -Wfloat-equal -Wno-unused-parameter)

add_executable(make_geo make_geo.cxx)
target_compile_features(make_geo PUBLIC cxx_std_23)
target_link_libraries(make_geo PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static)
target_compile_options(make_geo PRIVATE ${WARNINGS} $<$<CONFIG:Debug>:-g -O0> $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -fno-plt>
)
