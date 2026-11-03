# Still to reconsider

Every test, reproducible and script below predates tests/ and is up for retirement: replaced by a
simpler test in tests/, moved to examples/, or deleted. None of them is a reason not to rename or
change code; a file that breaks is fixed or retired. Strike an entry when it is settled, and name
its successor if there is one.

## Replaced (successor in tests/, oracle check pending)

- reproducibles_python/diffusion_convergence_elliptic.py -> tests/diffusion_elliptic.cxx

## Open

### reproducibles_python/

- reproducibles_python/advection_convergence_leVeque_hypergraph.py
- reproducibles_python/advection_convergence_parabolic.py
- reproducibles_python/advection_convergence_rotation.py
- reproducibles_python/advection_diffusion_reaction_convergence_elliptic.py
- reproducibles_python/advection_diffusion_reaction_hypergraph_convergence_elliptic.py
- reproducibles_python/advection_hypergraph.py
- reproducibles_python/bilaplacian_convergence_eigenvalue_approx.py
- reproducibles_python/bilaplacian_convergence_eigenvalue_comb.py
- reproducibles_python/bilaplacian_convergence_eigenvalue_newt.py
- reproducibles_python/bilaplacian_convergence_eigenvalue_shifted_inverse.py
- reproducibles_python/bilaplacian_convergence_elliptic.py
- reproducibles_python/bilaplacian_convergence_parabolic.py
- reproducibles_python/bilaplacian_hypergraph_convergence_eigenvalue_approx.py
- reproducibles_python/bilaplacian_hypergraph_convergence_eigenvalue_comb.py
- reproducibles_python/bilaplacian_hypergraph_convergence_eigenvalue_newt.py
- reproducibles_python/bilaplacian_hypergraph_convergence_eigenvalue_shifted_inverse.py
- reproducibles_python/bilaplacian_hypergraph_convergence_elliptic.py
- reproducibles_python/diffusion_convergence_eigenvalue_approx.py
- reproducibles_python/diffusion_convergence_eigenvalue_comb.py
- reproducibles_python/diffusion_convergence_eigenvalue_cross.py
- reproducibles_python/diffusion_convergence_eigenvalue_newt.py
- reproducibles_python/diffusion_convergence_eigenvalue_shifted_inverse.py
- reproducibles_python/diffusion_convergence_eigenvalue_shifted_inverse_cross.py
- reproducibles_python/diffusion_convergence_elliptic_ip.py
- reproducibles_python/diffusion_convergence_elliptic_postprocess.py
- reproducibles_python/diffusion_convergence_hypergraph_elliptic.py
- reproducibles_python/diffusion_convergence_hypergraph_elliptic_high_order.py
- reproducibles_python/diffusion_convergence_hypergraph_elliptic_high_order_postprocess.py
- reproducibles_python/diffusion_convergence_parabolic.py
- reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_approx.py
- reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_comb.py
- reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_newt.py
- reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_newt_cross.py
- reproducibles_python/diffusion_hypergraph_convergence_eigenvalue_shifted_inverse.py
- reproducibles_python/diffusion_hypergraph_convergence_elliptic.py
- reproducibles_python/diffusion_hypergraph_convergence_elliptic_ip.py
- reproducibles_python/diffusion_hypergraph_convergence_elliptic_postprocess.py
- reproducibles_python/movie_advection_hypergraph.py
- reproducibles_python/parameters/advection.hxx
- reproducibles_python/parameters/bilaplacian.hxx
- reproducibles_python/parameters/diffusion.hxx
- reproducibles_python/parameters/diffusion_advection_reaction.hxx
- reproducibles_python/parameters/timoshenko.hxx
- reproducibles_python/timoshenko_convergence_cross.py

### tests_c++/

- tests_c++/CMakeLists.txt
- tests_c++/diffusion_dimensions.cxx
- tests_c++/diffusion_floatings.cxx
- tests_c++/diffusion_plotting.cxx
- tests_c++/domain_hdf5.cxx
- tests_c++/plot_refined.cxx
- tests_c++/plot_simple.cxx
- tests_c++/plot_vtkhdf.cxx
- tests_c++/plot_vtkhdf_paraview.py
- tests_c++/point_functions.cxx
- tests_c++/qr_handedness.cxx
- tests_c++/stiffness_matrix.cxx

### tests_python/

- tests_python/CMakeLists.txt
- tests_python/diffusion_uniform.py
- tests_python/local_solvers/run.py
- tests_python/nanobind_bindings.py
- tests_python/operator_matrices.py
- tests_python/stiffness_matrix.py

### python/

- python/README.md
- python/check_vs_cython.py
- python/convergence.sh
- python/diffusion.py
- python/geo2h5.py
- python/hyperhdg.py
- python/operator_loop.py
- python/operator_loop.sh
- python/operator_loop_arrays.py
- python/operator_loop_check.py
- python/operator_loop_convergence.py

### examples/

- examples/CMakeLists.txt
- examples/diffusion_elliptic.cxx
- examples/parameters/diffusion.hxx

### cython/

- cython/elliptic.cfg
- cython/elliptic.pxd
- cython/elliptic.pyx
- cython/mass_approx_eigenvalue.cfg
- cython/mass_approx_eigenvalue.pxd
- cython/mass_approx_eigenvalue.pyx
- cython/nonlinear_eigenvalue.cfg
- cython/nonlinear_eigenvalue.pxd
- cython/nonlinear_eigenvalue.pyx
- cython/parabolic.cfg
- cython/parabolic.pxd
- cython/parabolic.pyx
- cython/shifted_inverse_eigenvalue.cfg
- cython/shifted_inverse_eigenvalue.pxd
- cython/shifted_inverse_eigenvalue.pyx

### import/

- import/HyperHDG/__init__.py
- import/HyperHDG/import_cxx/cmake.py
- import/HyperHDG/import_cxx/compile_prep.py
- import/HyperHDG/import_cxx/config.py
- import/HyperHDG/import_cxx/cpp.py
- import/HyperHDG/import_cxx/include.py
- import/HyperHDG/import_cxx/names.py
- import/HyperHDG/import_cxx/paths.py
- import/HyperHDG/scipy_compat.py

