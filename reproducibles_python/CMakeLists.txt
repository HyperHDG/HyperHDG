# The reproducibles run on precompiled nanobind modules: modules/*.cxx, one module per family of
# instantiations (global loop, local solver, topology), built here into
# ${CMAKE_BINARY_DIR}/reproducibles_python/ where common.py finds them (HYPERHDG_BUILD_DIR, or
# build/ next to the sources). One ctest per migrated script, in quick mode.

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
if(NOT TARGET nanobind-static)
  add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/nanobind.git ${CMAKE_BINARY_DIR}/nanobind)
endif()

file(GLOB module_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/modules/*.cxx)
foreach(source ${module_sources})
  get_filename_component(module ${source} NAME_WE)
  # target repro_<module> (examples/ has executables of the same names), file <module>.<abi>.so
  # NOMINSIZE: optimize for speed (-O3), not size (nanobind's default for release builds).
  nanobind_add_module(repro_${module} NOMINSIZE ${source})
  target_link_libraries(repro_${module} PRIVATE HyperHDG)
  target_include_directories(repro_${module} PRIVATE ${PROJECT_SOURCE_DIR})
  set_target_properties(repro_${module} PROPERTIES
    OUTPUT_NAME ${module} LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  list(APPEND modules repro_${module})
endforeach()
add_custom_target(reproducibles DEPENDS ${modules})

# Every script, in quick mode, checked against its expected records.
set(migrated
  diffusion_convergence_elliptic
  diffusion_convergence_elliptic_ip
  diffusion_convergence_elliptic_postprocess
  bilaplacian_convergence_elliptic
  advection_diffusion_reaction_convergence_elliptic
  diffusion_hypergraph_convergence_elliptic
  diffusion_hypergraph_convergence_elliptic_ip
  diffusion_hypergraph_convergence_elliptic_postprocess
  bilaplacian_hypergraph_convergence_elliptic
  advection_diffusion_reaction_hypergraph_convergence_elliptic
  diffusion_convergence_hypergraph_elliptic
  diffusion_convergence_hypergraph_elliptic_high_order
  diffusion_convergence_hypergraph_elliptic_high_order_postprocess
  diffusion_convergence_parabolic
  bilaplacian_convergence_parabolic
  advection_convergence_parabolic
  advection_convergence_rotation
  diffusion_convergence_eigenvalue_approx
  diffusion_convergence_eigenvalue_shifted_inverse
  diffusion_convergence_eigenvalue_newt
  diffusion_convergence_eigenvalue_comb
  diffusion_hypergraph_convergence_eigenvalue_approx
  diffusion_hypergraph_convergence_eigenvalue_shifted_inverse
  diffusion_hypergraph_convergence_eigenvalue_newt
  diffusion_hypergraph_convergence_eigenvalue_comb
  bilaplacian_convergence_eigenvalue_approx
  bilaplacian_convergence_eigenvalue_shifted_inverse
  bilaplacian_convergence_eigenvalue_newt
  bilaplacian_convergence_eigenvalue_comb
  bilaplacian_hypergraph_convergence_eigenvalue_approx
  bilaplacian_hypergraph_convergence_eigenvalue_shifted_inverse
  bilaplacian_hypergraph_convergence_eigenvalue_newt
  bilaplacian_hypergraph_convergence_eigenvalue_comb
  diffusion_convergence_eigenvalue_shifted_inverse_cross
  diffusion_hypergraph_convergence_eigenvalue_newt_cross
  diffusion_convergence_eigenvalue_cross
  advection_convergence_leVeque_hypergraph
  movie_advection_hypergraph
  timoshenko_convergence_cross
)
foreach(script ${migrated})
  add_test(
    NAME reproducible_${script}
    COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/${script}.py --quick --check
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
  # one BLAS thread: the dense local solves are small, threaded BLAS only spins
  set_tests_properties(reproducible_${script} PROPERTIES
    ENVIRONMENT "HYPERHDG_BUILD_DIR=${CMAKE_BINARY_DIR};OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
endforeach()
