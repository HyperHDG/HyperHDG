## Command-line tools. geo2h5 converts .geo domain files to .geo.h5 and needs HDF5.
if(HYPERHDG_HDF5_FOUND)
  add_executable(geo2h5 geo2h5.cxx)
  target_link_libraries(geo2h5 PRIVATE HyperHDG)

  # Every domain file of the repository converts (the files are written to the build tree).
  if(BUILD_TESTING)
    file(GLOB domain_files ${PROJECT_SOURCE_DIR}/domains/*.geo)
    foreach(domain_file ${domain_files})
      get_filename_component(domain ${domain_file} NAME_WE)
      add_test(NAME geo2h5_${domain}
               COMMAND geo2h5 ${domain_file} ${CMAKE_CURRENT_BINARY_DIR}/${domain}.geo.h5)
    endforeach()
  endif()
endif()
