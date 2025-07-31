if(NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/tensor_product_chain_complex.git/CMakeLists.txt)
  message(FATAL_ERROR "Git submodules missing: Either clone with --recursive or run 'git submodule update --init'")
endif()

add_subdirectory(tensor_product_chain_complex.git)
add_subdirectory(tensor_product_polynomials.git)
add_subdirectory(nanoflann)
add_subdirectory(fmtlog)
add_subdirectory(bxzstr)
add_subdirectory(KaHIP)
add_library(KaHIP INTERFACE)
target_include_directories(KaHIP INTERFACE
  KaHIP/lib/
  KaHIP/lib/io
  KaHIP/lib/partition
  KaHIP/lib/tools
)
add_subdirectory(CLI11)
add_subdirectory(METIS)
