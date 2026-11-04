# RECORDS: the test checks its record lines against expected/<name>.yaml (tests/records.hxx);
# HYPERHDG_RECORD=1 ctest -R <name> writes that file instead.
function(add_hyperhdg_test name source)
  cmake_parse_arguments(PARSE_ARGV 2 arg "RECORDS" "" "")
  add_executable(${name} ${source})
  # Debug builds keep the assertions but optimize: at -O0 the solves dominate the CI job.
  target_compile_options(${name} PRIVATE $<$<CONFIG:Debug>:-O2>)
  # The local matrices are SmallMats on the stack, and setting up one local operator of Diffusion
  # in 3D, degree 3, takes about 2.5 MB of it (A alone is 256 x 256 doubles), more than the 2 MB of
  # MinGW's default; 8 MB as on Linux. TODO: optimize the stack usage, by filling a heap workspace
  # of the loop in place instead of returning blocks and operators by value.
  target_link_options(${name} PRIVATE $<$<BOOL:${MINGW}>:-Wl,--stack,8388608>)
  target_link_libraries(${name} PRIVATE HyperHDG)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR})
  if(arg_RECORDS)
    target_compile_definitions(
      ${name} PRIVATE RECORDS_FILE="${CMAKE_CURRENT_SOURCE_DIR}/expected/${name}.yaml")
  endif()
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}")
endfunction()

# One executable per instantiation (dimension, degree): the build parallelizes over them, and ccache
# recompiles only the instantiations whose headers changed.
function(add_hyperhdg_ladder_test source)
  get_filename_component(test ${source} NAME_WE)
  foreach(dim 1 2 3)
    foreach(degree 1 2 3)
      add_hyperhdg_test(${test}_d${dim}_p${degree} ${source})
      target_compile_definitions(${test}_d${dim}_p${degree} PRIVATE DIM=${dim} DEGREE=${degree})
    endforeach()
  endforeach()
endfunction()

add_hyperhdg_test(cg cg.cxx)
add_hyperhdg_test(dense_la dense_la.cxx)
add_hyperhdg_test(generic generic.cxx)
add_hyperhdg_test(records records.cxx)
add_hyperhdg_ladder_test(diffusion_elliptic.cxx)

# Hyperedges of dimension 1 embedded in spaces of dimension 1 to 5 (domains/simplex_1_*.geo).
foreach(space_dim 1 2 3 4 5)
  add_hyperhdg_test(diffusion_embedded_s${space_dim} diffusion_embedded.cxx)
  target_compile_definitions(diffusion_embedded_s${space_dim} PRIVATE SPACE_DIM=${space_dim})
  set_tests_properties(diffusion_embedded_s${space_dim} PROPERTIES WORKING_DIRECTORY
                                                                   ${PROJECT_SOURCE_DIR})
endforeach()

# Single precision: one case, 2D degree 1, whose errors stay well above the round-off floor of the
# solve in float (degree 2 reaches it at n = 8: its error there depends on the rounding).
add_hyperhdg_test(diffusion_elliptic_float_d2_p1 diffusion_elliptic.cxx)
target_compile_definitions(diffusion_elliptic_float_d2_p1 PRIVATE DIM=2 DEGREE=1 REAL=float)

# DiffusionTensor (sum factorization) on the same convergence ladder.
foreach(dim 1 2 3)
  foreach(degree 1 2 3)
    add_hyperhdg_test(diffusion_tensor_d${dim}_p${degree} diffusion_elliptic.cxx)
    target_compile_definitions(diffusion_tensor_d${dim}_p${degree}
                               PRIVATE DIM=${dim} DEGREE=${degree} SOLVER=DiffusionTensor)
  endforeach()
endforeach()

# The Python side (nanobind is set up in the top-level CMakeLists.txt):
# - bindings.py: the NumPy bindings of Generic on a module built here (bindings_module.cxx).
# - hyperhdg_jit.py: python/hyperhdg.py on trivial modules, which it compiles just in time; a CMake
#   configure and build per module, so ctest -j starts it first and runs the others meanwhile.
if(NOT NOPYTHONTESTS AND Python_FOUND)
  nanobind_add_module(hyperhdg_bindings_test NOMINSIZE bindings_module.cxx)
  target_link_libraries(hyperhdg_bindings_test PRIVATE HyperHDG)
  target_include_directories(hyperhdg_bindings_test PRIVATE ${PROJECT_SOURCE_DIR})
  set_target_properties(hyperhdg_bindings_test PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                          ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME bindings_py COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/bindings.py)
  set_tests_properties(bindings_py PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR};${HYPERHDG_TEST_ENVIRONMENT}")

  add_test(NAME hyperhdg_jit_py
           COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/hyperhdg_jit.py
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(hyperhdg_jit_py PROPERTIES ENVIRONMENT "${HYPERHDG_TEST_ENVIRONMENT}"
                                                  COST 100)
endif()

# Plotting and domain files: their own runner (.github/workflows/plotting.yml), off by default.
option(HYPERHDG_PLOTTING_TESTS "Plotting and .geo.h5 tests (HDF5; pvbatch for the ParaView check)"
       OFF)
if(HYPERHDG_PLOTTING_TESTS)
  add_subdirectory(plotting)
endif()

# Quadruple precision end to end (std::float128_t, GCC >= 13). Off by default until the pinned TPP
# computes in the target precision (tensor_product_polynomials, branch quadrature-precision: the
# Gauss--Legendre rules, the shape functions at the points and the sums of the integrator); with
# the older TPP the error stops at about 1e-16 and the test fails.
option(HYPERHDG_QUAD_TESTS "Quadruple-precision test (std::float128_t; needs a current TPP)" OFF)
if(HYPERHDG_QUAD_TESTS)
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_FLAGS -std=c++23)
  check_cxx_source_compiles("#include <stdfloat>
    int main() { std::float128_t x = 1; return x > 0 ? 0 : 1; }" HYPERHDG_HAVE_FLOAT128)
  unset(CMAKE_REQUIRED_FLAGS)
  if(HYPERHDG_HAVE_FLOAT128)
    add_hyperhdg_test(quad quad.cxx)
    target_compile_features(quad PRIVATE cxx_std_23)
  else()
    message(WARNING "HYPERHDG_QUAD_TESTS: the compiler has no std::float128_t, quad left out.")
  endif()
endif()

# The Timoshenko beam solver on Generic (cross of beams, domains/cross.geo; run in the source tree),
# one executable per degree.
foreach(degree 1 2 3)
  add_hyperhdg_test(timoshenko_p${degree} timoshenko.cxx)
  target_compile_definitions(timoshenko_p${degree} PRIVATE DEGREE=${degree})
  set_tests_properties(timoshenko_p${degree} PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endforeach()
