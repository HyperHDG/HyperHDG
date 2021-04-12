enable_testing()


file (GLOB cpp_sources *.cxx)

foreach(ccfile ${cpp_sources})
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} ${LAPACK_LIBRARIES})
  add_test(NAME "${file}_cxx" COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()


if(NOT NOPYTHONTESTS)
  file (GLOB py_sources *.py)

  foreach(pyfile ${py_sources})
    get_filename_component(file ${pyfile} NAME_WE)
    add_test(NAME "${file}_py"
      COMMAND Python3::Interpreter ${pyfile}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples)
  endforeach()
endif()
