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

# The scripts already migrated to the precompiled modules; the others still go through the cython
# path and the weekly reproducibles workflow.
set(migrated
  diffusion_convergence_elliptic
)
foreach(script ${migrated})
  add_test(
    NAME reproducible_${script}
    COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/${script}.py --quick
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
  set_tests_properties(reproducible_${script} PROPERTIES
    ENVIRONMENT "HYPERHDG_BUILD_DIR=${CMAKE_BINARY_DIR}")
endforeach()
