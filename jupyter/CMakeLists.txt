enable_testing()

execute_process(COMMAND jupyter nbconvert --to script ${PROJECT_SOURCE_DIR}/jupyter/*.ipynb)

file (GLOB sources *.py)

foreach(pyfile ${sources})
  get_filename_component(file ${pyfile} NAME_WE)
  add_test(
    NAME ${file}_ipynb
    COMMAND Python3::Interpreter ${pyfile}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/jupyter
  )
endforeach()
