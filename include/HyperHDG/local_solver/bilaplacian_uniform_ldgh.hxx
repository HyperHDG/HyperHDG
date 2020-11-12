#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/quadrature_tensorial.hxx>
#include <HyperHDG/shape_fun_1d.hxx>
#include <HyperHDG/tensorial_shape_fun.hxx>
#include <algorithm>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief   Local solver for bilaplacian equation on uniform hypergraph.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t = double>
class BilaplacianUniform
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct empty_class
  {
  } data_type;

  // -----------------------------------------------------------------------------------------------
  // Public, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Dimension of hyper edge type that this object solves on.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   *
   * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
   * cooperating with this object.
   *
   * \retval  n_dofs        Number of global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node()
  {
    return 2 * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 2; }

 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = Hypercube<hyEdge_dimT>::pow(poly_deg + 1);
  /*!***********************************************************************************************
   * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  /*!***********************************************************************************************
   * \brief   Number of (local) degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = 2 * (hyEdge_dimT + 1) * n_shape_fct_;
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int system_dim = system_dimension();
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dim = node_system_dimension();
  /*!***********************************************************************************************
   * \brief   Assemble local matrix for the local solver.
   *
   * The local solver neither depends on the geometry, nor on global functions. Thus, its local
   * matrix is the same for all hyperedges and can be assembled once in the constructor. This is
   * done in this function.
   *
   * The function is static, since it is used in the constructor's initializer list.
   *
   * \param   tau           Penalty parameter for HDG.
   * \retval  loc_mat       Matrix of the local solver.
   ************************************************************************************************/
  static SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix(const lSol_float_t tau);
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;
  /*!***********************************************************************************************
   * \brief   Local matrix for the local solver.
   ************************************************************************************************/
  const SmallSquareMat<n_loc_dofs_, lSol_float_t> loc_mat_;
  /*!***********************************************************************************************
   * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
   ************************************************************************************************/
  const IntegratorTensorial<poly_deg, quad_deg, Gaussian, Legendre, lSol_float_t> integrator;

  // -----------------------------------------------------------------------------------------------
  // Private, internal functions for the local solver
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief  Assemble local right hand for the local solver.
   *
   * The right hand side needs the values of the global degrees of freedom. Thus, it needs to be
   * constructed individually for every hyperedge.
   *
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \retval  loc_rhs       Local right hand side of the locasl solver.
   ************************************************************************************************/
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs(
    const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values)
    const;
  /*!***********************************************************************************************
   * \brief   Solve local problem.
   *
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \retval  loc_sol       Solution of the local problem.
   ************************************************************************************************/
  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem(
    const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values)
    const
  {
    try
    {
      return (assemble_rhs(lambda_values) / loc_mat_).data();
    }
    catch (Wrapper::LAPACKexception& exc)
    {
      hy_assert(0 == 1, exc.what() << std::endl
                                   << "This can happen if quadrature is too inaccurate!");
      throw exc;
    }
  }
  /*!***********************************************************************************************
   * \brief   Evaluate primal variable at boundary.
   *
   * Function to evaluate primal variable of the solution. This function is needed to calculate
   * the local numerical fluxes.
   *
   * \param   coeffs        Coefficients of the local solution.
   * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
   ************************************************************************************************/
  inline std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> primal_at_boundary(
    const std::array<lSol_float_t, n_loc_dofs_>& coeffs) const;
  /*!***********************************************************************************************
   * \brief   Evaluate dual variable at boundary.
   *
   * Function to evaluate dual variable of the solution. This function is needed to calculate the
   * local numerical fluxes.
   *
   * \param   coeffs        Coefficients of the local solution.
   * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
   ************************************************************************************************/
  inline std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> dual_at_boundary(
    const std::array<lSol_float_t, 2 * (hyEdge_dimT + 1) * n_shape_fct_>& coeffs) const;

 public:
  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef lSol_float_t constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  BilaplacianUniform(const constructor_value_type& tau = 1.)
  : tau_(tau), loc_mat_(assemble_loc_matrix(tau))
  {
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
   * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
   * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
   * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
   *
   * \param   lambda_values Local part of vector x.
   * \param   time          Time --- this parameter is redundant for this local solver.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>
  numerical_flux_from_lambda(
    const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    std::array<lSol_float_t, n_loc_dofs_> coeffs = solve_local_problem(lambda_values);

    std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values,
      primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
        bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];

    return bdr_values;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector residual.
   *
   * \param   lambda_values Local part of vector x.
   * \param   time          Time --- this parameter is redundant for this local solver.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> numerical_flux_total(
    const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    std::array<lSol_float_t, n_loc_dofs_> coeffs = solve_local_problem(lambda_values);

    std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values,
      primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
        bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];

    return bdr_values;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local squared L2 error.
   *
   * \param   lambda_values The values of the skeletal variable's coefficients.
   * \param   time          Time --- this parameter is redundant for this local solver.
   * \retval  vec_b         Local part of vector b.
   ************************************************************************************************/
  lSol_float_t calc_L2_error_squared(
    const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  function_values   Function values at quadrature points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dim> bulk_values(
    const std::array<abscissa_float_t, sizeT>& abscissas,
    const input_array_t& lambda_values,
    const lSol_float_t time = 0.) const;
  /*!***********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   boundary_number   number of the boundary on which to evaluate the function.
   * \retval  function_values   Function values at quadrature points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>, node_system_dim>
  lambda_values(const std::array<abscissa_float_t, sizeT>& abscissas,
                const input_array_t& lambda_values,
                const unsigned int boundary_number) const;

};  // end of class BilaplacianUniform

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF BilaplacianUniform
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
SmallSquareMat<BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_loc_dofs_,
               lSol_float_t>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::assemble_loc_matrix(
  const lSol_float_t tau)
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  const IntegratorTensorial<poly_deg, quad_deg, Gaussian, Legendre, lSol_float_t> integrator;
  lSol_float_t integral;

  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dimT>(i, j);
      local_mat(hyEdge_dimT * n_shape_fct_ + i, n_dofs_lap + hyEdge_dimT * n_shape_fct_ + j) -=
        integral;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) += integral;
        local_mat(n_dofs_lap + dim * n_shape_fct_ + i, n_dofs_lap + dim * n_shape_fct_ + j) +=
          integral;
      }

      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        // Integral_element - nabla phi_i \vec phi_j dx
        // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
        integral = integrator.template integrate_vol_Dphiphi<hyEdge_dimT>(i, j, dim);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) -= integral;
        local_mat(dim * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) -= integral;
        local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i,
                  n_dofs_lap + dim * n_shape_fct_ + j) -= integral;
        local_mat(n_dofs_lap + dim * n_shape_fct_ + i,
                  n_dofs_lap + hyEdge_dimT * n_shape_fct_ + j) -= integral;

        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 1);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) += integral;
        local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i,
                  n_dofs_lap + dim * n_shape_fct_ + j) +=
          integral;  // Removing to enforce Neumann zero!
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
        local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i,
                  n_dofs_lap + hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 0);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) -= integral;
        local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i,
                  n_dofs_lap + dim * n_shape_fct_ + j) -=
          integral;  // Removing to enforce Neumann zero!
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
        local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i,
                  n_dofs_lap + hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
      }
    }
  }

  return local_mat;
}  // end of BernoulliBendingBeam::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline SmallVec<BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_loc_dofs_,
                lSol_float_t>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::assemble_rhs(
  const std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values)
  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  lSol_float_t integral;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;

  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == 2 * n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2 * dim + 0][j] * integral;
        right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * lambda_values[2 * dim + 0][j] * integral;
        right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i] +=
          lambda_values[2 * dim + 0][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * lambda_values[2 * dim + 0][n_shape_bdr_ + j] * integral;

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        right_hand_side[dim * n_shape_fct_ + i] -= lambda_values[2 * dim + 1][j] * integral;
        right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * lambda_values[2 * dim + 1][j] * integral;
        right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i] -=
          lambda_values[2 * dim + 1][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * lambda_values[2 * dim + 1][n_shape_bdr_ + j] * integral;
      }
    }
  }

  return right_hand_side;
}  // end of BernoulliBendingBeam::assemble_rhs

// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline std::array<
  std::array<lSol_float_t,
             2 * BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::primal_at_boundary(
  const std::array<lSol_float_t, n_loc_dofs_>& coeffs) const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)
    bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        bdr_values[2 * dim + 0][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
        bdr_values[2 * dim + 0][n_shape_bdr_ + j] +=
          coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] * integral;

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2 * dim + 1][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
        bdr_values[2 * dim + 1][n_shape_bdr_ + j] +=
          coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] * integral;
      }
    }
  }

  return bdr_values;
}  // end of BilaplacianUniform::primal_at_boundary

// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline std::array<
  std::array<lSol_float_t,
             2 * BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::dual_at_boundary(
  const std::array<lSol_float_t, 2 * (hyEdge_dimT + 1) * n_shape_fct_>& coeffs) const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array<std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)
    bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        bdr_values[2 * dim + 0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
        bdr_values[2 * dim + 0][n_shape_bdr_ + j] -=
          coeffs[n_dofs_lap + dim * n_shape_fct_ + i] * integral;

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2 * dim + 1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
        bdr_values[2 * dim + 1][n_shape_bdr_ + j] +=
          coeffs[n_dofs_lap + dim * n_shape_fct_ + i] * integral;
      }
    }
  }

  return bdr_values;
}  // end of BilaplacianUniform::dual_at_boundary

// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
           BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::system_dim>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::bulk_values(
  const std::array<abscissa_float_t, sizeT>& abscissas,
  const input_array_t& lambda_values,
  const lSol_float_t time) const
{
  std::array<lSol_float_t, n_loc_dofs_> coefficients = solve_local_problem(lambda_values);
  std::array<lSol_float_t, n_loc_dofs_ / 2> first_half_of_coefficients;
  std::copy_n(coefficients.begin(), n_loc_dofs_ / 2, first_half_of_coefficients.begin());
  TensorialShapeFunctionEvaluation<hyEdge_dimT, lSol_float_t, Legendre, poly_deg, sizeT,
                                   abscissa_float_t>
    evaluation(abscissas);
  return evaluation
    .template evaluate_linear_combination_in_all_tensorial_points<system_dimension()>(
      first_half_of_coefficients);

}  // end of BilaplacianUniform::bulk_values

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
           BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::node_system_dim>
BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::lambda_values(
  const std::array<abscissa_float_t, sizeT>& abscissas,
  const input_array_t& lambda_values,
  const unsigned int boundary_number) const
{
  return std::array<
    std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
    BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::node_system_dimension()>();
}

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

}  // namespace LocalSolver
