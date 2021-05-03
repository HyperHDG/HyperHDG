enable_testing()

file (GLOB sources *.py)

foreach(pyfile ${sources})
  get_filename_component(file ${pyfile} NAME_WE)
  add_test(
    NAME ${file}_py
    COMMAND Python3::Interpreter ${pyfile}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
endforeach()
