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
# nanobind is set up in the top-level CMakeLists.txt.
if(NOT NOPYTHONTESTS AND Python_FOUND)
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
  # A CMake configure and build: started first by ctest -j, so the others run meanwhile.
  set_tests_properties(diffusion_elliptic_nb_jit_py PROPERTIES
    ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}" COST 100)
  # A network of beams (domains/grid_7x7.geo), solved and plotted; also compiled just in time.
  add_test(NAME beam_network_nb_jit_py
           COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/beam_network_nb_jit.py
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(beam_network_nb_jit_py PROPERTIES
    ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}" COST 100)
endif()
