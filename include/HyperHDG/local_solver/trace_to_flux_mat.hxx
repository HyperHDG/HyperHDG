#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <array>

/*!*************************************************************************************************
 * \brief   Define member function \c trace_to_flux_mat of a local solver.
 *
 * The defined function returns the local condensed matrix, i.e., the matrix representation of
 * \c trace_to_flux: Column \c face * n_glob_dofs_per_node() + \c dof is the result of
 * \c trace_to_flux applied to the unit vector that is one for degree of freedom \c dof of face
 * \c face (and its rows are ordered analogously). The local matrix is assembled and LU factorized
 * only once, and the factorization is reused for all unit vectors. The results are the same as the
 * ones of calling \c trace_to_flux for all unit vectors.
 *
 * The macro is to be used within the class definition of a local solver that provides
 * - \c hyEdge_dimT, \c lSol_float_t, \c n_loc_dofs_, and \c n_glob_dofs_per_node(),
 * - \c assemble_rhs_from_lambda(lambda_values, hyper_edge) returning a \c SmallVec,
 * - \c flux_from_coeffs(lambda_values_in, coeffs, lambda_values_out, hyper_edge), which evaluates
 *   \c trace_to_flux for given local solution \c coeffs (the solution of the local problem),
 * and whose \c trace_to_flux solves the local problem with the matrix \c loc_mat_expr (which may
 * use \c hyper_edge and \c time).
 *
 * The local condensed matrix may only depend on the hyperedge's geometry and node descriptor, the
 * local solver's parameters, and \c time, but not on the data stored at the hyperedge
 * (\c hyper_edge.data), just like \c trace_to_flux.
 *
 * \c assemble_rhs_from_lambda and \c flux_from_coeffs are called for every local unit vector, and
 * usually evaluate boundary integrals in every call. For local solvers whose local matrix is
 * cheap compared to these integrals, a hand-written \c trace_to_flux_mat that evaluates the
 * integrals only once can be considerably faster (cf. \c Diffusion::trace_to_flux_mat).
 **************************************************************************************************/
#define HYPERHDG_LOCAL_TRACE_TO_FLUX_MAT(loc_mat_expr)                                            \
  template <typename hyEdgeT>                                                                     \
  SmallSquareMat<2 * hyEdge_dimT * n_glob_dofs_per_node(), lSol_float_t> trace_to_flux_mat(       \
    hyEdgeT& hyper_edge, const lSol_float_t time = 0.) const                                      \
  {                                                                                               \
    constexpr unsigned int n_dofs = n_glob_dofs_per_node();                                       \
    SmallSquareMat<2 * hyEdge_dimT * n_dofs, lSol_float_t> loc_mat;                               \
    std::array<std::array<lSol_float_t, n_dofs>, 2 * hyEdge_dimT> lambda_in, lambda_out;          \
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)                                   \
      lambda_in[face].fill(0.);                                                                   \
                                                                                                  \
    SmallSquareMat<n_loc_dofs_, lSol_float_t> lu_mat = loc_mat_expr;                              \
    std::array<int, n_loc_dofs_> ipiv;                                                            \
    try                                                                                           \
    {                                                                                             \
      Wrapper::lapack_factorize<n_loc_dofs_, lSol_float_t>(lu_mat.data(), ipiv);                  \
    }                                                                                             \
    catch (Wrapper::LAPACKexception & exc)                                                        \
    {                                                                                             \
      hy_assert(0 == 1, exc.what() << std::endl                                                   \
                                   << "This can happen if quadrature is too inaccurate!");        \
      throw exc;                                                                                  \
    }                                                                                             \
                                                                                                  \
    for (unsigned int face_j = 0; face_j < 2 * hyEdge_dimT; ++face_j)                             \
      for (unsigned int dof_j = 0; dof_j < n_dofs; ++dof_j)                                       \
      {                                                                                           \
        lambda_in[face_j][dof_j] = 1.;                                                            \
        SmallVec<n_loc_dofs_, lSol_float_t> coeffs =                                              \
          assemble_rhs_from_lambda(lambda_in, hyper_edge);                                        \
        Wrapper::lapack_solve_factored<n_loc_dofs_, 1, lSol_float_t>(lu_mat.data(), ipiv,         \
                                                                     coeffs.data());              \
        flux_from_coeffs(lambda_in, coeffs, lambda_out, hyper_edge);                              \
        for (unsigned int face_i = 0; face_i < 2 * hyEdge_dimT; ++face_i)                         \
          for (unsigned int dof_i = 0; dof_i < n_dofs; ++dof_i)                                   \
            loc_mat(face_i* n_dofs + dof_i, face_j * n_dofs + dof_j) = lambda_out[face_i][dof_i]; \
        lambda_in[face_j][dof_j] = 0.;                                                            \
      }                                                                                           \
                                                                                                  \
    return loc_mat;                                                                               \
  }                                                                                               \
  static_assert(true, "")
