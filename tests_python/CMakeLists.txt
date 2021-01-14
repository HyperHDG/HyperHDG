enable_testing()

file (GLOB sources *.py)

foreach(pyfile ${sources})
  get_filename_component(file ${pyfile} NAME_WE)
  add_test(NAME ${file}
    COMMAND Python3::Interpreter -B ${pyfile}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests_python
  )
endforeach()