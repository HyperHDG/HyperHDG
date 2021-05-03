enable_testing()

file (GLOB sources *.ipynb)

foreach(ipynbfile ${sources})
  get_filename_component(file ${ipynbfile} NAME_WE)
  add_test(
    NAME ${file}_ipynb
    COMMAND
      jupyter nbconvert --execute --to script ${ipynbfile}
      Python3::Interpreter ${file}.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/jupyter
  )
endforeach()
