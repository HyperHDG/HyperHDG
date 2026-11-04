function(add_hyperhdg_test ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC HyperHDG)
  add_test(NAME ${file} COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endfunction()

add_hyperhdg_test(diffusion_dimensions.cxx)
add_hyperhdg_test(diffusion_floatings.cxx)
