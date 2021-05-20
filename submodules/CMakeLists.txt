if(NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/tensor_product_chain_complex.git/CMakeLists.txt)
  message(FATAL_ERROR "Git submodules missing: Either clone with --recursive or run 'git submodule update --init'")
endif()

add_subdirectory(tensor_product_chain_complex.git)
add_subdirectory(tensor_product_polynomials.git)
