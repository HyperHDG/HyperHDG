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

# Compile and run every local solver once (cf. local_solvers/run.py).
foreach(case
    Diffusion_Elliptic Diffusion_Elliptic_hypergraph DiffusionIP_Elliptic
    DiffusionPostprocess_Elliptic DiffusionAdvectionReaction_Elliptic Bilaplacian_Elliptic
    BilaplacianUniform_Elliptic TimoshenkoBeam_Elliptic DiffusionParab_Parabolic
    BilaplacianParab_Parabolic AdvectionParab_Parabolic DiffusionEigs_NonlinearEigenvalue
    BilaplacianEigs_NonlinearEigenvalue DiffusionEigs_ShiftedInverseEigenvalue
    BilaplacianEigs_ShiftedInverseEigenvalue Diffusion_MassApproxEigenvalue
    Bilaplacian_MassApproxEigenvalue)
  add_test(
    NAME local_solver_${case}_py
    COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/local_solvers/run.py ${case}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
endforeach()
