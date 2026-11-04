# Plotting (VTU, VTKHDF) and domain files (.geo.h5), label "plotting": run by
# .github/workflows/plotting.yml when the plotting or file-reading code changes, not by the main
# runner. The VTKHDF and .geo.h5 tests need HDF5, the ParaView check pvbatch (-DPVBATCH_EXECUTABLE
# or on PATH). The tests write into output/ of the source tree; the target plotting_tests builds
# them all.
add_custom_target(plotting_tests)

function(add_plotting_test name)
  add_executable(${name} ${name}.cxx)
  target_compile_options(${name} PRIVATE $<$<CONFIG:Debug>:-O2>)
  target_link_libraries(${name} PRIVATE HyperHDG)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR})
  add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(${name} PROPERTIES LABELS plotting
                                          ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
  add_dependencies(plotting_tests ${name})
endfunction()

add_plotting_test(plot_simple)
add_plotting_test(plot_refined)
add_plotting_test(diffusion_plotting)
add_plotting_test(plot_exact)
set(plot_files "")
foreach(dims 1-1 1-2 1-3 2-2 2-3 3-3)
  list(APPEND plot_files ${PROJECT_SOURCE_DIR}/output/plot_exact-${dims}.vtu)
endforeach()
set(plot_writers plot_exact)

if(HYPERHDG_HDF5_FOUND)
  add_plotting_test(domain_hdf5)
  add_plotting_test(plot_vtkhdf)
  foreach(dims 1-1 1-2 1-3 2-2 2-3 3-3)
    list(APPEND plot_files ${PROJECT_SOURCE_DIR}/output/plot_vtkhdf-${dims}.vtkhdf)
  endforeach()
  list(APPEND plot_writers plot_vtkhdf)
else()
  message(WARNING "HDF5 not found: the VTKHDF and .geo.h5 tests are left out.")
endif()

# ParaView reads every plot back (values checked at every point) and renders it to
# output/plots/<name>.png.
find_program(PVBATCH_EXECUTABLE pvbatch)
# Kitware's headless (OSMesa) builds need --mesa: their launcher then puts lib/mesa, where
# libOSMesa finds libglapi, on LD_LIBRARY_PATH. A pvbatch of a distribution takes no options.
set(PVBATCH_OPTIONS "" CACHE STRING "Options of pvbatch, e.g. --mesa for Kitware's OSMesa build")
if(PVBATCH_EXECUTABLE)
  add_test(NAME plot_paraview
           COMMAND ${PVBATCH_EXECUTABLE} ${PVBATCH_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/check_plots.py
                   ${plot_files}
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(plot_paraview PROPERTIES LABELS plotting DEPENDS "${plot_writers}")
else()
  message(WARNING "pvbatch not found: the ParaView check of the plots is left out.")
endif()
