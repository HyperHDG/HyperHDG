# Every example is also a ctest, so that the examples are known to work.

function(add_hyperdg_example ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC HyperHDG)
  target_include_directories(${file} PUBLIC ${PROJECT_SOURCE_DIR})
  add_test(NAME "${file}_cxx" COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties("${file}_cxx" PROPERTIES ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
endfunction()

add_hyperdg_example(diffusion_elliptic.cxx)

# The nanobind examples: a module built here (diffusion_elliptic_nb_module.cxx) and the script
# that imports it, and the script that compiles its module just in time (python/hyperhdg.py).
if(NOT NOPYTHONTESTS)
  # nanobind looks up Python as "Python"; use the interpreter the top level found as "Python3".
  if(Python3_EXECUTABLE AND NOT Python_EXECUTABLE)
    set(Python_EXECUTABLE ${Python3_EXECUTABLE})
  endif()
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module)
  if(Python_FOUND)
    if(NOT TARGET nanobind-static)
      add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/nanobind.git ${CMAKE_BINARY_DIR}/nanobind)
    endif()
    nanobind_add_module(diffusion_elliptic_nb NOMINSIZE diffusion_elliptic_nb_module.cxx)
    target_link_libraries(diffusion_elliptic_nb PRIVATE HyperHDG)
    target_include_directories(diffusion_elliptic_nb PRIVATE ${PROJECT_SOURCE_DIR})
    set_target_properties(diffusion_elliptic_nb PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                           ${CMAKE_CURRENT_BINARY_DIR})

    add_test(NAME diffusion_elliptic_nb_precompiled_py
             COMMAND Python::Interpreter
                     ${CMAKE_CURRENT_SOURCE_DIR}/diffusion_elliptic_nb_precompiled.py
             WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    set_tests_properties(diffusion_elliptic_nb_precompiled_py PROPERTIES
      ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR};${HYPERHDG_TEST_ENVIRONMENT}")
    add_test(NAME diffusion_elliptic_nb_jit_py
             COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/diffusion_elliptic_nb_jit.py
             WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    set_tests_properties(diffusion_elliptic_nb_jit_py PROPERTIES
      ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
  endif()
endif()
