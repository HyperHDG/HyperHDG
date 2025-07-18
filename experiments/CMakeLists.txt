add_executable(make_geo make_geo.cxx)
target_compile_features(make_geo PUBLIC cxx_std_23)
target_link_libraries(make_geo PRIVATE nanoflann::nanoflann PRIVATE fmtlog-static)
