# Benchmarks with nanobench (vendored in extern/), not part of CI: see README.md. Configure with
# -DCMAKE_BUILD_TYPE=Release -DHYPERHDG_PERFORMANCE_TESTS=ON; run with ctest -L performance or the
# executables directly. Every benchmark is serial (one BLAS/OpenMP thread, RUN_SERIAL), its
# results go to stdout and to output/performance/<name>.json.
add_library(hyperhdg_nanobench STATIC nanobench.cxx)
target_include_directories(hyperhdg_nanobench PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern)
add_custom_target(performance)

function(add_hyperhdg_benchmark name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE HyperHDG hyperhdg_nanobench)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR})
  add_test(NAME performance_${name} COMMAND ${name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(performance_${name} PROPERTIES LABELS performance RUN_SERIAL TRUE
                                                      ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
  add_dependencies(performance ${name})
endfunction()

# One executable per dimension: the instantiations of a dimension compile in parallel to the other.
foreach(dim 2 3)
  add_hyperhdg_benchmark(local_operators_d${dim} local_operators.cxx)
  target_compile_definitions(local_operators_d${dim} PRIVATE DIM=${dim})
endforeach()
add_hyperhdg_benchmark(parameters parameters.cxx)
if(TARGET SuiteSparse::CHOLMOD)
  add_hyperhdg_benchmark(linear_solvers linear_solvers.cxx)
  target_link_libraries(linear_solvers PRIVATE SuiteSparse::CHOLMOD)
else()
  message(WARNING "CHOLMOD (SuiteSparse) not found: the benchmark linear_solvers is left out.")
endif()
