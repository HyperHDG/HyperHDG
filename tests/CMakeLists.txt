function(add_hyperhdg_test name source)
  add_executable(${name} ${source})
  # Debug builds keep the assertions but optimize: at -O0 the solves dominate the CI job.
  target_compile_options(${name} PRIVATE $<$<CONFIG:Debug>:-O2>)
  target_link_libraries(${name} PRIVATE HyperHDG)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR})
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
endfunction()

# One executable per instantiation (dimension, degree): the build parallelizes over them, and ccache
# recompiles only the instantiations whose headers changed.
function(add_hyperhdg_ladder_test source)
  get_filename_component(test ${source} NAME_WE)
  foreach(dim 1 2 3)
    foreach(degree 1 2 3)
      add_hyperhdg_test(${test}_d${dim}_p${degree} ${source})
      target_compile_definitions(${test}_d${dim}_p${degree} PRIVATE DIM=${dim} DEGREE=${degree})
    endforeach()
  endforeach()
endfunction()

add_hyperhdg_test(cg cg.cxx)
add_hyperhdg_test(dense_la dense_la.cxx)
add_hyperhdg_test(generic generic.cxx)
add_hyperhdg_ladder_test(diffusion_elliptic.cxx)

# Hyperedges of dimension 1 embedded in spaces of dimension 1 to 5 (domains/simplex_1_*.geo).
foreach(space_dim 1 2 3 4 5)
  add_hyperhdg_test(diffusion_embedded_s${space_dim} diffusion_embedded.cxx)
  target_compile_definitions(diffusion_embedded_s${space_dim} PRIVATE SPACE_DIM=${space_dim})
  set_tests_properties(diffusion_embedded_s${space_dim} PROPERTIES WORKING_DIRECTORY
                                                                   ${PROJECT_SOURCE_DIR})
endforeach()

# Single precision: one case, 2D degree 2, whose errors stay well above float round-off.
add_hyperhdg_test(diffusion_elliptic_float_d2_p2 diffusion_elliptic.cxx)
target_compile_definitions(diffusion_elliptic_float_d2_p2 PRIVATE DIM=2 DEGREE=2 REAL=float)

# The Python side (nanobind is set up in the top-level CMakeLists.txt):
# - bindings.py: the NumPy bindings of Generic on a module built here (bindings_module.cxx).
# - hyperhdg_jit.py: python/hyperhdg.py on trivial modules, which it compiles just in time; a CMake
#   configure and build per module, so ctest -j starts it first and runs the others meanwhile.
if(NOT NOPYTHONTESTS AND Python_FOUND)
  nanobind_add_module(hyperhdg_bindings_test NOMINSIZE bindings_module.cxx)
  target_link_libraries(hyperhdg_bindings_test PRIVATE HyperHDG)
  target_include_directories(hyperhdg_bindings_test PRIVATE ${PROJECT_SOURCE_DIR})
  set_target_properties(hyperhdg_bindings_test PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                          ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME bindings_py COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/bindings.py)
  set_tests_properties(bindings_py PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR};${HYPERHDG_TEST_ENVIRONMENT}")

  add_test(NAME hyperhdg_jit_py
           COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/hyperhdg_jit.py
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(hyperhdg_jit_py PROPERTIES ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}"
                                                  COST 100)
endif()

# Plotting and domain files: their own runner (.github/workflows/plotting.yml), off by default.
option(HYPERHDG_PLOTTING_TESTS "Plotting and .geo.h5 tests (HDF5; pvbatch for the ParaView check)"
       OFF)
if(HYPERHDG_PLOTTING_TESTS)
  add_subdirectory(plotting)
endif()
