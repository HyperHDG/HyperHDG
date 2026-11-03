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
add_hyperhdg_test(generic generic.cxx)
add_hyperhdg_ladder_test(diffusion_elliptic.cxx)
