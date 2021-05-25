function(add_hyperdg_example ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC HyperHDG)
  target_include_directories(${file} PUBLIC ${PROJECT_SOURCE_DIR})
  add_test(NAME "${file}_cxx" COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endfunction()

add_hyperdg_example(diffusion_elliptic.cxx)


if(NOT NOPYTHONTESTS)
  file (GLOB py_sources *.py)

  foreach(pyfile ${py_sources})
    get_filename_component(file ${pyfile} NAME_WE)
    add_test(NAME "${file}_py"
      COMMAND Python3::Interpreter ${pyfile}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endforeach()
endif()
