function(add_hyperhdg_test ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC HyperHDG)
  add_test(NAME ${file} COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endfunction()

add_hyperhdg_test(diffusion_dimensions.cxx)
add_hyperhdg_test(diffusion_floatings.cxx)
add_hyperhdg_test(diffusion_plotting.cxx)
add_hyperhdg_test(plot_refined.cxx)
add_hyperhdg_test(plot_simple.cxx)
add_hyperhdg_test(point_functions.cxx)
add_hyperhdg_test(qr_handedness.cxx)
add_hyperhdg_test(stiffness_matrix.cxx)

# Using factorize_once for non-uniform local matrices must be detected (and abort) in debug mode.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_test(NAME stiffness_matrix_misuse COMMAND stiffness_matrix misuse
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(stiffness_matrix_misuse PROPERTIES WILL_FAIL TRUE)
endif()
