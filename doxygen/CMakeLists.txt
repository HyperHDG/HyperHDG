# The documentation as .github/workflows/make_doxygen.yml builds it (doxygen/Doxyfile; TPP and TPCC
# first, so that their tag files exist), output in <build>/doxygen/html.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(doxygen_env ${CMAKE_COMMAND} -E env DOXYGEN_OUTPUT_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(
    doxygen
    # Doxygen creates one directory level only, not html/tpp without html
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html
    COMMAND ${doxygen_env} ${DOXYGEN_EXECUTABLE} doxygen/Doxyfile.tpp
    COMMAND ${doxygen_env} ${DOXYGEN_EXECUTABLE} doxygen/Doxyfile.tpcc
    COMMAND ${doxygen_env} ${DOXYGEN_EXECUTABLE} doxygen/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Building the Doxygen documentation"
    VERBATIM
  )
endif()
