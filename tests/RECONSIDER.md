# Still to reconsider

Every test, reproducible and script below predates tests/ and is up for retirement: replaced by a
simpler test in tests/, moved to examples/, or deleted. None of them is a reason not to rename or
change code; a file that breaks is fixed or retired. Strike an entry when it is settled, and name
its successor if there is one.

## Replaced (successor in tests/, oracle check pending)

- reproducibles_python/diffusion_convergence_elliptic.py -> tests/diffusion_elliptic.cxx
- tests_c++/point_functions.cxx, tests_c++/qr_handedness.cxx -> tests/dense_la.cxx
- tests_c++/diffusion_dimensions.cxx, its Diffusion part -> tests/diffusion_embedded.cxx (the
  DiffusionUniform part stays legacy; its do_test_standard expects 0.5 at the free nodes, the
  solution of -u'' = 0 with Kirchhoff's law is 1 / (1 + sqrt 2), which Elliptic and Generic both
  give with Dirichlet data from the parameters)

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
- tests_c++/diffusion_floatings.cxx (DiffusionUniform in float, stays legacy with it; float on
  Generic + Diffusion: tests/diffusion_elliptic_float_d2_p2)
- tests_c++/diffusion_plotting.cxx
- tests_c++/domain_hdf5.cxx
- tests_c++/plot_refined.cxx
- tests_c++/plot_simple.cxx
- tests_c++/plot_vtkhdf.cxx
- tests_c++/plot_vtkhdf_paraview.py

### tests_python/

- tests_python/CMakeLists.txt
- tests_python/local_solvers/run.py
- tests_python/nanobind_bindings.py

### python/

- python/README.md
- python/geo2h5.py
- python/hyperhdg.py

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

