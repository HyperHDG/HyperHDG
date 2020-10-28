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
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo    Implement local solvers with template parameters in_vector and out_vector.
 *
 * This class contains the local solver for Poisson's equation, i.e.,
 * \f[
 *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t = double>
class DiffusionUniform
{
 public:
  typedef struct empty_class
  {
  } data_type;

  typedef lSol_float_t solver_float_t;

  // -----------------------------------------------------------------------------------------------
  // Public, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief Dimension of hyper edge type that this object solves on.
   *
   * \todo  Why is this not just called dimension?
   *        -> E.g. in elasticity there are two important dimensions, the one of the hyperedge and
   *        the one of the space. Thus, elasticity can return both dimensions, while this class
   *        only returns the relevant hyperedge dimension.
   * \todo  The original brief referred to the internal variable only. It should be the other way
   *        round: this function is the main access to this number.
   *        -> I agree, this was not on purpose and I have to check for this in other classes!
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   *
   * \todo  Why are these called global degrees of freedom and not just `n_dofs_per_node()`?
   *        -> In Elasticity, there are two types of dofs per node. The one that come from outside
   *        (they are space_dim - dimensional) and the ones that are relevant for the local
   *        problem (and therefore hyEdge_dimT - dimensional). Thus, there is a discrimination
   *        between global and local amount per dofs in local solvers.
   *
   *
   * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
   * cooperating with this object.
   *
   * \retval  n_dofs        Number of global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node()
  {
    return Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_value_dimension() { return 1U; }
  /*!***********************************************************************************************
   * \brief Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }

  static constexpr unsigned int node_system_dimension() { return 1; }

 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1);
  /*!***********************************************************************************************
   * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node();
  /*!***********************************************************************************************
   * \brief   Number of (local) degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = (hyEdge_dimT + 1) * n_shape_fct_;

  static constexpr unsigned int system_dim = system_dimension();

  static constexpr unsigned int node_system_dim = node_system_dimension();

  /*!***********************************************************************************************
   * \brief  Assemble local matrix for the local solver.
   *
   * The local solver neither depends on the geometry, nor on global functions. Thus, its local
   * matrix is the same for all hyperedges and can be assembled once in the constructor. This is
   * done in this function.
   *
   * The function is static inline, since it is used in the constructor's initializer list.
   *
   * \param   tau           Penalty parameter for HDG.
   * \retval  loc_mat       Matrix of the local solver.
   ************************************************************************************************/
  static SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix(const lSol_float_t tau);
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;
  /*!*********************************************************************************************
   * \brief   Local matrix for the local solver.
   **********************************************************************************************/
  const SmallSquareMat<n_loc_dofs_, lSol_float_t> loc_mat_;

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
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values) const;

  /*!***********************************************************************************************
   * \brief  Solve local problem.
   *
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \retval  loc_sol       Solution of the local problem.
   ************************************************************************************************/
  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem(
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values) const
  {
    try
    {
      return (assemble_rhs(lambda_values) / loc_mat_).data();
    }
    catch (LAPACKexception& exc)
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
  inline std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> primal_at_boundary(
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
  inline std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> dual_at_boundary(
    const std::array<lSol_float_t, (hyEdge_dimT + 1) * n_shape_fct_>& coeffs) const;

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
  DiffusionUniform(const constructor_value_type& tau = 1.)
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
   * \param   time          Time.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> numerical_flux_from_lambda(
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    std::array<lSol_float_t, n_loc_dofs_> coeffs = solve_local_problem(lambda_values);

    std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values,
      primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
        bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];

    return bdr_values;
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
   * \param   time          Time.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> numerical_flux_total(
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    return numerical_flux_from_lambda(lambda_values, time);
  }

  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \param   lambda_values The values of the skeletal variable's coefficients.
   * \param   time          Time.
   * \retval  vec_b         Local part of vector b.
   ************************************************************************************************/
  lSol_float_t calc_L2_error_squared(
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t time = 0.) const
  {
    return 0.;
  }

  template <typename abscissa_float_t, std::size_t abscissas_sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>, system_dim>
  bulk_values(const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
              const input_array_t& lambda_values,
              const lSol_float_t time = 0.) const;

  /*!***********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary
   *
   * \tparam  absc_float_t  Floating type for the abscissa values.
   * \tparam  abscissas_sizeT         Size of the array of array of abscissas.
   * \param   abscissas     Abscissas of the supporting points.
   * \param   lambda_values The values of the skeletal variable's coefficients.
   * \param   boundary_number number of the boundary on which to evaluate the function.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t abscissas_sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(abscissas_sizeT)>,
             node_system_dim>

  lambda_values(const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
                const input_array_t& lambda_values,
                const unsigned int boundary_number) const;

};  // end of class DiffusionUniform

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF DiffusionUniform
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
SmallSquareMat<DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_loc_dofs_,
               lSol_float_t>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::assemble_loc_matrix(
  const lSol_float_t tau)
{
  const IntegratorTensorial<poly_deg, quad_deg, Gaussian, Legendre, lSol_float_t> integrator;
  lSol_float_t integral;

  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dimT>(i, j);
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) += integral;

      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        // Integral_element - nabla phi_i \vec phi_j dx
        // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
        integral = integrator.template integrate_vol_Dphiphi<hyEdge_dimT>(i, j, dim);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) -= integral;
        local_mat(dim * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) -= integral;

        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 1);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) += integral;
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 0);
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) -= integral;
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) += tau * integral;
      }
    }
  }

  return local_mat;
}  // end of DiffusionUniform::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline SmallVec<DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_loc_dofs_,
                lSol_float_t>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::assemble_rhs(
  const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values) const
{
  lSol_float_t integral;

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;

  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == n_shape_bdr_,
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

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        right_hand_side[dim * n_shape_fct_ + i] -= lambda_values[2 * dim + 1][j] * integral;
        right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * lambda_values[2 * dim + 1][j] * integral;
      }
    }
  }

  return right_hand_side;
}  // end of DiffusionUniform::assemble_rhs

// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline std::array<
  std::array<lSol_float_t,
             DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::primal_at_boundary(
  const std::array<lSol_float_t, n_loc_dofs_>& coeffs) const
{
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;
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

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2 * dim + 1][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
      }
    }
  }

  return bdr_values;
}  // end of DiffusionUniform::primal_at_boundary

// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
inline std::array<
  std::array<lSol_float_t,
             DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::dual_at_boundary(
  const std::array<lSol_float_t, (hyEdge_dimT + 1) * n_shape_fct_>& coeffs) const
{
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;
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

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2 * dim + 1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
      }
    }
  }

  return bdr_values;
}  // end of DiffusionUniform::dual_at_boundary

// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
template <typename abscissa_float_t, std::size_t abscissas_sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
           DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::system_dim>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::bulk_values(
  const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
  const input_array_t& lambda_values,
  const lSol_float_t time) const
{
  std::array<lSol_float_t, n_loc_dofs_> coefficients = solve_local_problem(lambda_values);

  TensorialShapeFunctionEvaluation<hyEdge_dimT, lSol_float_t, Legendre, poly_deg, abscissas_sizeT,
                                   abscissa_float_t>
    evaluation(abscissas);
  return evaluation
    .template evaluate_linear_combination_in_all_tensorial_points<system_dimension()>(coefficients);
}  // end of DiffusionUniform::bulk_values

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t>
template <typename abscissa_float_t, std::size_t abscissas_sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(abscissas_sizeT)>,
           DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::node_system_dim>
DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t>::lambda_values(
  const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
  const input_array_t& lambda_values,
  const unsigned int boundary_number) const
{
  TensorialShapeFunctionEvaluation<hyEdge_dimT - 1, lSol_float_t, Legendre, poly_deg,
                                   abscissas_sizeT, abscissa_float_t>
    evaluation(abscissas);
  return evaluation
    .template evaluate_linear_combination_in_all_tensorial_points<node_system_dimension()>(
      lambda_values[boundary_number]);
}

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

}  // namespace LocalSolver
