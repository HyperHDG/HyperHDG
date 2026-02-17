if(NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/tensor_product_chain_complex.git/CMakeLists.txt)
  message(FATAL_ERROR "Git submodules missing: Either clone with --recursive or run 'git submodule update --init'")
endif()

add_subdirectory(tensor_product_chain_complex.git)
add_subdirectory(tensor_product_polynomials.git)
add_subdirectory(nanoflann)
# HACK: overwrite default setting of NOMPI
set(NOMPI ON CACHE BOOL "disable all targets that depend on MPI (kaffpaE, ParHIP)")
add_subdirectory(KaHIP)

# HACK: ignore -Wunused in KaHIP
get_directory_property(KAHIP_TARGETS DIRECTORY KaHIP BUILDSYSTEM_TARGETS)
if (NOT NOMPI AND PARHIP)
  get_directory_property(PARHIP_TARGETS DIRECTORY KaHIP/parallel/parallel_src BUILDSYSTEM_TARGETS)
  get_directory_property(MODKAHIP_TARGETS DIRECTORY KaHIP/parallel/modified_kahip BUILDSYSTEM_TARGETS)
endif()
foreach(target ${KAHIP_TARGETS} ${PARHIP_TARGETS} ${MODKAHIP_TARGETS})
  target_compile_options(${target} PRIVATE -Wno-unused)
endforeach()
