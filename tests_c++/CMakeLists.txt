enable_testing()

file (GLOB sources *.cxx)

foreach(ccfile ${sources})
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} ${LAPACK_LIBRARIES})
  add_test(NAME ${file} COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
