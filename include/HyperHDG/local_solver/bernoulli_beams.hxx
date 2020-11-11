#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/quadrature_tensorial.hxx>
#include <HyperHDG/shape_fun_1d.hxx>
#include <HyperHDG/tensorial_shape_fun.hxx>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief   Local solver for the equation that governs the lengthening of an elastic beam.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t = double,
          typename diffusion_sol_t =
            DiffusionUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t> >
class LengtheningBeam
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct empty_class
  {
  } data_type;
  /*!***********************************************************************************************
   *  \brief  Define floating type the local solver uses for use of external classses / functions.
   ************************************************************************************************/
  typedef lSol_float_t solver_float_t;
  /*!***********************************************************************************************
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
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
    return 2 * space_dim * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return space_dim * 2; }

 private:
  /*!***********************************************************************************************
   * \brief   The diffusion solver that solves the locally defined PDE.
   ************************************************************************************************/
  const diffusion_sol_t diffusion;
  /*!***********************************************************************************************
   * \brief   Do the pretprocessing to transfer global to local dofs.
   ************************************************************************************************/
  template <class hyEdgeT>
  inline std::array<std::array<double, diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  node_dof_to_edge_dof(
    const std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>& lambda,
    hyEdgeT& hyper_edge) const
  {
    std::array<std::array<double, diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result;
    hy_assert(result.size() == 2, "Only implemented in one dimension!");
    for (unsigned int i = 0; i < result.size(); ++i)
    {
      hy_assert(result[i].size() == 1, "Only implemented in one dimension!");
      result[i].fill(0.);
    }

    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(1);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        result[i][0] += normal_vector[dim] * lambda[i][dim];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Do the postprocessing to transfer local to global dofs.
   ************************************************************************************************/
  template <class hyEdgeT>
  inline std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  edge_dof_to_node_dof(const std::array<std::array<double, diffusion_sol_t::n_glob_dofs_per_node()>,
                                        2 * hyEdge_dimT>& lambda,
                       hyEdgeT& hyper_edge) const
  {
    hy_assert(diffusion_sol_t::n_glob_dofs_per_node() == 1, "This should be 1!");
    std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result;
    for (unsigned int i = 0; i < result.size(); ++i)
      result[i].fill(0.);
    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(1);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        result[i][dim] = normal_vector[dim] * lambda[i][0];

    return result;
  }

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
  LengtheningBeam(const constructor_value_type& tau = 1.) : diffusion(tau) {}
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    HyperEdge that is considered.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_from_lambda(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                              2 * hyEdge_dimT>& lambda_values,
                             hyEdgeT& hyper_edge,
                             const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "Elastic graphs must be graphs, not hypergraphs!");
    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge);

    if constexpr (not_uses_geometry<
                    diffusion_sol_t,
                    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                               2 * hyEdge_dimT>(
                      std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                                 2 * hyEdge_dimT>&)>::value)
      lambda = diffusion.numerical_flux_from_lambda(lambda);
    else
      lambda = diffusion.numerical_flux_from_lambda(lambda, hyper_edge);

    return edge_dof_to_node_dof(lambda, hyper_edge);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to residual.
   *
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    HyperEdge that is considered.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vecAx         Local part of vector A * x - b.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_total(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                        2 * hyEdge_dimT>& lambda_values,
                       hyEdgeT& hyper_edge,
                       const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "Elastic graphs must be graphs, not hypergraphs!");
    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge);

    if constexpr (not_uses_geometry<
                    diffusion_sol_t,
                    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                               2 * hyEdge_dimT>(
                      std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                                 2 * hyEdge_dimT>&)>::value)
      lambda = diffusion.numerical_flux_total(lambda);
    else
      lambda = diffusion.numerical_flux_total(lambda, hyper_edge);

    return edge_dof_to_node_dof(lambda, hyper_edge);
  }
  /*!***********************************************************************************************
   * \brief   Local squared contribution to the L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  vec_b             Local part of vector b.
   ************************************************************************************************/
  template <class hyEdgeT>
  lSol_float_t calc_L2_error_squared(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    lSol_float_t error = 0.;
    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge);

    if constexpr (not_uses_geometry<
                    diffusion_sol_t,
                    std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                               2 * hyEdge_dimT>(
                      std::array<std::array<lSol_float_t, diffusion_sol_t::n_glob_dofs_per_node()>,
                                 2 * hyEdge_dimT>&)>::value)
      error = diffusion.calc_L2_error_squared(lambda);
    else
      error = diffusion.calc_L2_error_squared(lambda, hyper_edge);

    return error;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time.
   * \retval  func_values       Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<
    std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
    LengtheningBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, lSol_float_t>::system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    std::array<
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      LengtheningBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, lSol_float_t>::system_dimension()>
      result;

    auto bulk = diffusion.bulk_values(abscissas, node_dof_to_edge_dof(lambda_values, hyper_edge));
    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(1);

    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] = bulk[1][q] * normal_vector[dim];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary.
   *
   * \tparam  absc_float_t    Floating type for the abscissa values.
   * \tparam  abscissas_sizeT Size of the array of array of abscissas.
   * \param   abscissas       Abscissas of the supporting points.
   * \param   lambda_values   The values of the skeletal variable's coefficients.
   * \param   boundary_number Number of the boundary on which to evaluate the function.
   * \retval  func_values     Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
             node_system_dimension()>
  lambda_values(const std::array<abscissa_float_t, sizeT>& abscissas,
                const input_array_t& lambda_values,
                const unsigned int boundary_number) const
  {
    TensorialShapeFunctionEvaluation<hyEdge_dimT - 1, lSol_float_t, Legendre, poly_deg, sizeT,
                                     abscissa_float_t>
      evaluation(abscissas);
    return evaluation
      .template evaluate_linear_combination_in_all_tensorial_points<node_system_dimension()>(
        lambda_values[boundary_number]);
  }
};  // end of class LengtheningBeam

/*!*************************************************************************************************
 * \brief   Local solver for the equation that governs the bending of an elastic Bernoulli beam.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t = double,
          typename bilaplacian_sol_t =
            BilaplacianUniform<hyEdge_dimT, poly_deg, quad_deg, lSol_float_t> >
class BernoulliBendingBeam
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct empty_class
  {
  } data_type;
  /*!***********************************************************************************************
   *  \brief  Define floating type the local solver uses for use of external classses / functions.
   ************************************************************************************************/
  typedef lSol_float_t solver_float_t;
  /*!***********************************************************************************************
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
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
    return 2 * space_dim * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return space_dim * 2; }

 private:
  /*!***********************************************************************************************
   * \brief   The bilaplacian solver that solves the locally defined PDE.
   ************************************************************************************************/
  const bilaplacian_sol_t bilaplacian_solver;
  /*!***********************************************************************************************
   * \brief   Do the preprocessing to transfer global to local dofs.
   ************************************************************************************************/
  template <class hyEdgeT>
  inline std::array<std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  node_dof_to_edge_dof(
    const std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>& lambda,
    hyEdgeT& hyper_edge,
    const unsigned int outer_index) const
  {
    std::array<std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      result;
    hy_assert(result.size() == 2, "Only implemented in one dimension!");
    for (unsigned int i = 0; i < result.size(); ++i)
    {
      hy_assert(result[i].size() == 2, "Only implemented in one dimension!");
      result[i].fill(0.);
    }

    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(outer_index);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        result[i][0] += normal_vector[dim] * lambda[i][dim];
        result[i][1] += normal_vector[dim] * lambda[i][space_dim + dim];
      }

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Do the postprocessing to transfer local to global dofs.
   ************************************************************************************************/
  template <class hyEdgeT>
  inline std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  edge_dof_to_node_dof(
    const std::array<std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                     2 * hyEdge_dimT>& lambda,
    hyEdgeT& hyper_edge,
    const unsigned int outer_index) const
  {
    hy_assert(bilaplacian_sol_t::n_glob_dofs_per_node() == 2, "This should be 1*2!");
    std::array<std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result;
    for (unsigned int i = 0; i < result.size(); ++i)
      result[i].fill(0.);
    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(outer_index);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        result[i][dim] = normal_vector[dim] * lambda[i][0];
        result[i][space_dim + dim] = normal_vector[dim] * lambda[i][1];
      }

    return result;
  }

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
  BernoulliBendingBeam(const constructor_value_type& tau = 1.) : bilaplacian_solver(tau) {}
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    HyperEdge that is considered.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_from_lambda(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                              2 * hyEdge_dimT>& lambda_values,
                             hyEdgeT& hyper_edge,
                             const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "The beam must be one-dimensional!");
    std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda;

    std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result, aux;
    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      result[i].fill(0.);

    for (unsigned int dim = 0; dim < space_dim - hyEdge_dimT; ++dim)
    {
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge, dim);

      if constexpr (
        not_uses_geometry<
          bilaplacian_sol_t,
          std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                     2 * hyEdge_dimT>(
            std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                       2 * hyEdge_dimT>&)>::value)
        lambda = bilaplacian_solver.numerical_flux_from_lambda(lambda);
      else
        lambda = bilaplacian_solver.numerical_flux_from_lambda(lambda, hyper_edge);

      aux = edge_dof_to_node_dof(lambda, hyper_edge, dim);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
          result[i][j] += aux[i][j];
    }

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to residual.
   *
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    HyperEdge that is considered.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vecAx         Local part of vector A * x - b.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_total(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                        2 * hyEdge_dimT>& lambda_values,
                       hyEdgeT& hyper_edge,
                       const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "The beam must be one-dimensional!");
    std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda;

    std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result, aux;
    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      result[i].fill(0.);

    for (unsigned int dim = 0; dim < space_dim - hyEdge_dimT; ++dim)
    {
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge, dim);

      if constexpr (
        not_uses_geometry<
          bilaplacian_sol_t,
          std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                     2 * hyEdge_dimT>(
            std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                       2 * hyEdge_dimT>&)>::value)
        lambda = bilaplacian_solver.numerical_flux_total(lambda);
      else
        lambda = bilaplacian_solver.numerical_flux_total(lambda, hyper_edge);

      aux = edge_dof_to_node_dof(lambda, hyper_edge, dim);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
          result[i][j] += aux[i][j];
    }

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Local squared contribution to the L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  vec_b             Local part of vector b.
   ************************************************************************************************/
  template <class hyEdgeT>
  lSol_float_t calc_L2_error_squared(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    lSol_float_t error = 0.;
    std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
      lambda;
    for (unsigned int dim = 0; dim < space_dim - hyEdge_dimT; ++dim)
    {
      lambda = node_dof_to_edge_dof(lambda_values, hyper_edge, dim);

      if constexpr (
        not_uses_geometry<
          bilaplacian_sol_t,
          std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                     2 * hyEdge_dimT>(
            std::array<std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>,
                       2 * hyEdge_dimT>&)>::value)
        error += bilaplacian_solver.calc_L2_error_squared(lambda);
      else
        error += bilaplacian_solver.calc_L2_error_squared(lambda, hyper_edge);
    }

    return error;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time.
   * \retval  func_values       Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
             BernoulliBendingBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, lSol_float_t>::
               system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
      values;
    for (unsigned int i = 0; i < values.size(); ++i)
      values[i].fill(0.);

    for (unsigned int dim_on = 0; dim_on < space_dim - hyEdge_dimT; ++dim_on)
    {
      auto bulk = bilaplacian_solver.bulk_values(
        abscissas, node_dof_to_edge_dof(lambda_values, hyper_edge, dim_on));
      Point<space_dim, lSol_float_t> normal_vector =
        (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(dim_on);

      for (unsigned int dim = 0; dim < values.size(); ++dim)
        for (unsigned int q = 0; q < values[dim].size(); ++q)
          values[dim][q] += bulk[1][q] * normal_vector[dim];
    }

    return values;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary.
   *
   * \tparam  absc_float_t    Floating type for the abscissa values.
   * \tparam  abscissas_sizeT Size of the array of array of abscissas.
   * \param   abscissas       Abscissas of the supporting points.
   * \param   lambda_values   The values of the skeletal variable's coefficients.
   * \param   boundary_number Number of the boundary on which to evaluate the function.
   * \retval  func_values     Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
             node_system_dimension()>
  lambda_values(const std::array<abscissa_float_t, sizeT>& abscissas,
                const input_array_t& lambda_values,
                const unsigned int boundary_number) const
  {
    TensorialShapeFunctionEvaluation<hyEdge_dimT - 1, lSol_float_t, Legendre, poly_deg, sizeT,
                                     abscissa_float_t>
      evaluation(abscissas);
    return evaluation
      .template evaluate_linear_combination_in_all_tensorial_points<node_system_dimension()>(
        lambda_values[boundary_number]);
  }
};  // end of class BernoulliBendingBeam

/*!*************************************************************************************************
 * \brief   Local solver for the equation that governs the bending and change of length of an
 *          elastic Bernoulli beam.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          typename lSol_float_t = double>
class LengtheningBernoulliBendingBeam
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct empty_class
  {
  } data_type;
  /*!***********************************************************************************************
   *  \brief  Define floating type the local solver uses for use of external classses / functions.
   ************************************************************************************************/
  typedef lSol_float_t solver_float_t;
  /*!***********************************************************************************************
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
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
    return 2 * space_dim * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return space_dim * 2; }

 private:
  /*!***********************************************************************************************
   * \brief   The lengthening beam solver that does the lengthening of the beam.
   ************************************************************************************************/
  const LengtheningBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, lSol_float_t> len_beam;
  /*!***********************************************************************************************
   * \brief   The bending beam solver that does the bending of the beam.
   ************************************************************************************************/
  const BernoulliBendingBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, lSol_float_t> ben_beam;

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
  LengtheningBernoulliBendingBeam(const constructor_value_type& tau = 1.)
  : len_beam(tau), ben_beam(tau)
  {
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    Hyperedge.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vecAx         Local part of vector A * x.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_from_lambda(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                              2 * hyEdge_dimT>& lambda_values,
                             hyEdgeT& hyper_edge,
                             const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "A beam must be one-dimensional!");
    std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result, aux;

    result = len_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);
    aux = ben_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
        result[i][j] += aux[i][j];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to residual.
   *
   * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values Lambda values.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  vec_b         Local part of vector A * x - b.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>
  numerical_flux_total(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                        2 * hyEdge_dimT>& lambda_values,
                       hyEdgeT& hyper_edge,
                       const lSol_float_t time = 0.) const
  {
    static_assert(hyEdge_dimT == 1, "A beam must be one-dimensional!");
    std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT> result, aux;

    result = len_beam.numerical_flux_total(lambda_values, hyper_edge);
    aux = ben_beam.numerical_flux_total(lambda_values, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
        result[i][j] += aux[i][j];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Local squared contribution to the L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  vec_b             Local part of vector b.
   ************************************************************************************************/
  template <class hyEdgeT>
  lSol_float_t calc_L2_error_squared(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    lSol_float_t error = 0.;
    error += len_beam.calc_L2_error_squared(lambda_values, hyper_edge);
    error += ben_beam.calc_L2_error_squared(lambda_values, hyper_edge);
    return error;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time.
   * \retval  func_values       Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
      result, auxiliary;
    result = len_beam.bulk_values(abscissas, lambda_values, hyper_edge);
    auxiliary = ben_beam.bulk_values(abscissas, lambda_values, hyper_edge);

    for (unsigned int i = 0; i < system_dimension(); ++i)
      for (unsigned int j = 0; j < Hypercube<hyEdge_dimT>::pow(sizeT); ++j)
        result[i][j] += auxiliary[i][j];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary.
   *
   * \tparam  absc_float_t    Floating type for the abscissa values.
   * \tparam  abscissas_sizeT Size of the array of array of abscissas.
   * \param   abscissas       Abscissas of the supporting points.
   * \param   lambda_values   The values of the skeletal variable's coefficients.
   * \param   boundary_number Number of the boundary on which to evaluate the function.
   * \retval  func_values     Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
             node_system_dimension()>
  lambda_values(const std::array<abscissa_float_t, sizeT>& abscissas,
                const input_array_t& lambda_values,
                const unsigned int boundary_number) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
               node_system_dimension()>
      result = len_beam.lambda_values(abscissas, lambda_values, boundary_number);
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>,
               node_system_dimension()>
      auxiliary = ben_beam.lambda_values(abscissas, lambda_values, boundary_number);

    for (unsigned int i = 0; i < node_system_dimension(); ++i)
      for (unsigned int j = 0; j < Hypercube<hyEdge_dimT - 1>::pow(sizeT); ++j)
        result[i][j] += auxiliary[i][j];

    return result;
  }
};  // end of class LengtheningBernoulliBendingBeam

}  // namespace LocalSolver
