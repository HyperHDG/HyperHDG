# One executable per instantiation (dimension, degree): the build parallelizes over them, and ccache
# recompiles only the instantiations whose headers changed.
function(add_hyperhdg_ladder_test source)
  get_filename_component(test ${source} NAME_WE)
  foreach(dim 1 2 3)
    foreach(degree 1 2 3)
      set(name ${test}_d${dim}_p${degree})
      add_executable(${name} ${source})
      target_compile_definitions(${name} PRIVATE DIM=${dim} DEGREE=${degree})
      # Debug builds keep the assertions but optimize: at -O0 the solves dominate the CI job.
      target_compile_options(${name} PRIVATE $<$<CONFIG:Debug>:-O2>)
      target_link_libraries(${name} PRIVATE HyperHDG)
      target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR})
      add_test(NAME ${name} COMMAND ${name})
    endforeach()
  endforeach()
endfunction()

add_hyperhdg_ladder_test(diffusion_elliptic.cxx)
