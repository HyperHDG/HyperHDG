#pragma once  // Ensure that file is included only once in a single compilation.

// TODO: write the integrator function where the func can be vector valued in diff wave!! -> andreas

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <iostream>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <tuple>

template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave0
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};  // end of struct DiffusionParametersDefault

namespace LocalSolver
{

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
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
          template <unsigned int, typename> typename parametersT = TestTimoWave0,
          typename lSol_float_t = double>
class TimoshenkoWave
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of node elements, especially with respect to nodal shape functions.
   ************************************************************************************************/
  struct node_element
  {
    typedef std::tuple<TPP::ShapeFunction<
      TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1>>>
      functions;
  };
  /*!***********************************************************************************************
   *  \brief  Define how errors are evaluated.
   ************************************************************************************************/
  struct error_def
  {
    /*!*********************************************************************************************
     *  \brief  Define the typename returned by function errors.
     **********************************************************************************************/
    typedef std::array<lSol_float_t, 1U> error_t;
    /*!*********************************************************************************************
     *  \brief  Define how initial error is generated.
     **********************************************************************************************/
    static error_t initial_error()
    {
      std::array<lSol_float_t, 1U> summed_error;
      summed_error.fill(0.);
      return summed_error;
    }
    /*!*********************************************************************************************
     *  \brief  Define how local errors should be accumulated.
     **********************************************************************************************/
    static error_t sum_error(error_t& summed_error, const error_t& new_error)
    {
      for (unsigned int k = 0; k < summed_error.size(); ++k)
        summed_error[k] += new_error[k];
      return summed_error;
    }
    /*!*********************************************************************************************
     *  \brief  Define how global errors should be postprocessed.
     **********************************************************************************************/
    static error_t postprocess_error(error_t& summed_error)
    {
      for (unsigned int k = 0; k < summed_error.size(); ++k)
        summed_error[k] = std::sqrt(summed_error[k]);
      return summed_error;
    }
  };
  /*!***********************************************************************************************
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   *
   * This number must be equal to HyperNodeFactory::n_glob_dofs_per_node()() of the HyperNodeFactory
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
  static constexpr unsigned int system_dimension() { return 6*space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 6*space_dim; }

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
  static constexpr unsigned int n_loc_dofs_ = 6 * space_dim * n_shape_fct_;
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int system_dim = system_dimension();

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;
  const lSol_float_t theta_;
  const lSol_float_t delta_t_;

  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT>>,
    lSol_float_t>
    integrator;

  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT>>,
    lSol_float_t>
    vec_integrator;

  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
 public:
  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef std::vector<double> constructor_value_type;

  struct data_type
  {
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> u_old, v_old, r_old, s_old, flux_u, flux_r, n_old, m_old, flux_v, flux_s;
  };
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  // NOTE: tau, theta, delta_t
  TimoshenkoWave(const constructor_value_type& vals = std::vector(3, 1.)) : tau_(vals[0]),
    theta_(vals[1]), delta_t_(vals[2]) {}

  template <typename point_t, typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_vol_phivecfunccomp_beam(const unsigned int i,
                                                                         std::array<int, n_comps> comps, geom_t& geom, const lSol_float_t time) const
  {
    std::array<lSol_float_t, n_comps> ret;
    for (unsigned int j = 0; j < n_comps; j++) {
      ret[j] = integrator::template integrate_vol_phivecfunccomp<
        point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(i, comps[j], geom,
                                                                         time);
    }
    return ret;
  }


  template <typename point_t, typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_vol_phivecfunccomp_beam_avg(const unsigned int i,
                                                                         std::array<int, n_comps> comps, geom_t& geom, const lSol_float_t time) const
  {
    std::array<lSol_float_t, n_comps> ret;
    for (unsigned int j = 0; j < n_comps; j++) {
      ret[j] = theta_ * integrator::template integrate_vol_phivecfunccomp<
        point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(i, comps[j], geom,
                                                                         time);
      ret[j] += (1-theta_) * integrator::template integrate_vol_phivecfunccomp<
        point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(i, comps[j], geom,
                                                                         time-delta_t_);

    }
    return ret;
  }

  template <typename point_t,
            typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_bdr_phivecfunccomp_beam(const unsigned int i,
                                               const unsigned int bdr,
                                               std::array<int, n_comps> comps,
                                               geom_t& geom,
                                               const lSol_float_t time = 0.) const
  {
   std::array<lSol_float_t, n_comps> ret;
   for (unsigned int j = 0; j < n_comps; j++) {
     ret[j] = integrator::template integrate_bdr_phivecfunccomp<
       point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(
                 i, bdr, comps[j], geom, time);
   }

    return ret;
  }



  template <typename point_t,
            typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_bdr_phivecfunccomp_beam_avg(const unsigned int i,
                                               const unsigned int bdr,
                                               std::array<int, n_comps> comps,
                                               geom_t& geom,
                                               const lSol_float_t time = 0.) const
  {
   std::array<lSol_float_t, n_comps> ret;
   for (unsigned int j = 0; j < n_comps; j++) {
     ret[j] = theta_ * integrator::template integrate_bdr_phivecfunccomp<
       point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(
                 i, bdr, comps[j], geom, time);
     ret[j] += (1-theta_) * integrator::template integrate_bdr_phivecfunccomp<
       point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(
                 i, bdr, comps[j], geom, time-delta_t_);

   }

    return ret;
  }


  template <typename hyEdgeT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix(
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const;

  template <typename hyEdgeT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix_a(
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const;

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_lambda(
    const SmallMatT& lambda_values,
    hyEdgeT& hyper_edge) const;

  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs(
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const;

  template <class hyEdgeT, typename SmallMatT>
  inline SmallMatT glob_dof_to_loc_dof(
    const SmallMatT& glob_dof,
    hyEdgeT& hyper_edge) const
  {
    SmallMatT loc_dof;

    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        loc_dof[i] += normal_vec[dim] * glob_dof[i+dim*n_shape_fct_];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          loc_dof[i+(1+ind)*n_shape_fct_] += normal_vec[dim] * glob_dof[i+n_shape_fct_*dim];
        }
    }

    return loc_dof;
  }

  template <class hyEdgeT, typename SmallMatT>
  inline SmallMatT loc_dof_to_glob_dof(
    const SmallMatT& loc_dof,
    hyEdgeT& hyper_edge) const
  {
    SmallMatT glob_dof;
    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        glob_dof[i + dim * n_shape_fct_] += normal_vec[dim] * loc_dof[i];

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          glob_dof[i + dim * n_shape_fct_] += normal_vec[dim] * loc_dof[i + (1 + ind) * n_shape_fct_];
    }
    return glob_dof;
  }

  template <class hyEdgeT, typename SmallMatT>
  inline std::array<std::array<double, 2 * space_dim>, 2 * hyEdge_dimT> node_dof_to_edge_dof(
    const SmallMatT& glob_lambda,
    hyEdgeT& hyper_edge) const
  {
    std::array<std::array<double, 2 * space_dim>, 2 * hyEdge_dimT> loc_lambda;
    hy_assert(loc_lambda.size() == 2, "Only implemented in one dimension!");
    for (unsigned int i = 0; i < loc_lambda.size(); ++i)
    {
      hy_assert(loc_lambda[i].size() == 6, "Only implemented in one dimension!");
      loc_lambda[i].fill(0.);
    }

    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        loc_lambda[i][0] += normal_vec[dim] * glob_lambda[i][dim];
        loc_lambda[i][0 + space_dim] += normal_vec[dim] * glob_lambda[i][dim + space_dim];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          loc_lambda[i][1 + ind] += normal_vec[dim] * glob_lambda[i][dim];
          loc_lambda[i][1 + ind + space_dim] += normal_vec[dim] * glob_lambda[i][space_dim + dim];
        }
    }

    return loc_lambda;
  }

  template <class hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  inline SmallMatOutT& edge_dof_to_node_dof(const SmallMatInT& loc_lambda,
                                            SmallMatOutT& glob_lambda,
                                            hyEdgeT& hyper_edge) const
  {
    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        glob_lambda[i][dim] += normal_vec[dim] * loc_lambda[i][0];
        glob_lambda[i][dim + space_dim] += normal_vec[dim] * loc_lambda[i][space_dim];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          glob_lambda[i][dim] += normal_vec[dim] * loc_lambda[i][1 + ind];
          glob_lambda[i][dim + space_dim] += normal_vec[dim] * loc_lambda[i][1 + ind + space_dim];
        }
    }

    return glob_lambda;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> solve_local_problem(const SmallMatT& lambda_values,
                                                                 const unsigned int solution_type,
                                                                 hyEdgeT& hyper_edge,
                                                                 const lSol_float_t time) const
  {
    try
    {
      SmallVec<n_loc_dofs_, lSol_float_t> rhs;
      if (solution_type == 0)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
      else if (solution_type == 1)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge) +
              assemble_rhs_from_global_rhs(hyper_edge, time);
      else
        hy_assert(0 == 1, "This has not been implemented!");
      // std::cout << "-- solve_local" << std::endl;
      // std::cout << rhs << std::endl;
      return rhs / assemble_loc_matrix(hyper_edge, time);
    }
    catch (Wrapper::LAPACKexception& exc)
    {
      std::cout << hyper_edge.geometry.area() << std::endl;
      hy_assert(0 == 1, exc.what() << std::endl
                                   << "This can happen if quadrature is too inaccurate!");
      throw exc;
    }
  }

  template <typename hyEdgeT>
  inline SmallMat<2 * hyEdge_dimT, 2 * space_dim * n_shape_bdr_, lSol_float_t>
  extract_fluxes_from_coeffs(const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
                             hyEdgeT& hyper_edge) const
  {
    SmallMat<2 * hyEdge_dimT, 2 * space_dim * n_shape_bdr_, lSol_float_t> bdr_values;
    lSol_float_t integral;

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
            bdr_values(face, j + dim) +=
              hyper_edge.geometry.local_normal(face).operator[](0) * integral *
                coeffs[dim * n_shape_fct_ + i] +
              tau_ * coeffs[(2 * space_dim + dim) * n_shape_fct_ + i] * integral;
        }

    return bdr_values;
  }

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_to_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time = 0.) const
  {
    // std::cout << "------------- trace_to_flux" << std::endl;

    hy_assert(lambda_values_in.size() == lambda_values_out.size() &&
                lambda_values_in.size() == 2 * hyEdge_dimT,
              "Both matrices must be of same size which corresponds to the number of faces!");
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      hy_assert(
        lambda_values_in[i].size() == lambda_values_out[i].size() &&
          lambda_values_in[i].size() == n_glob_dofs_per_node(),
        "Both matrices must be of same size which corresponds to the number of dofs per face!");

    SmallMatInT lambda_in = lambda_values_in;

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i) {
      if (hyper_edge.node_descriptor[i] & (1<<6)) continue;  // static-only flag
      for (unsigned int j = 0; j < 2*space_dim; j++)
        if (hyper_edge.node_descriptor[i] & (1<<j))
          lambda_in[i][j] = 0.;
    }

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_in, hyper_edge);

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl;

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, 0U, hyper_edge, time);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = result(i, j) - tau_ * lambda_values_loc[i][j];

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << result(i,j) << " ";
    // std::cout << std::endl;

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl << std::endl;

    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i) {
      if (hyper_edge.node_descriptor[i] & (1<<6)) continue;  // static-only flag
      for (unsigned int j = 0; j < 2*space_dim; j++)
        if (hyper_edge.node_descriptor[i] & (1<<j))
          lambda_values_out[i][j] = 0.;
    }

    return lambda_values_out;
  }

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& residual_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time = 0.) const
  {
    // std::cout << "------------- residual_flux" << std::endl;

    hy_assert(lambda_values_in.size() == lambda_values_out.size() &&
                lambda_values_in.size() == 2 * hyEdge_dimT,
              "Both matrices must be of same size which corresponds to the number of faces!");
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      hy_assert(
        lambda_values_in[i].size() == lambda_values_out[i].size() &&
          lambda_values_in[i].size() == n_glob_dofs_per_node(),
        "Both matrices must be of same size which corresponds to the number of dofs per face!");

    SmallMatInT lambda_in = lambda_values_in;

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i) {
      if (hyper_edge.node_descriptor[i] & (1<<6)) continue;  // static-only flag
      for (unsigned int j = 0; j < 2*space_dim; j++)
        if (hyper_edge.node_descriptor[i] & (1<<j))
          lambda_in[i][j] = 0.;
    }

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_in, hyper_edge);

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, 1U, hyper_edge, time);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);
    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = result(i, j) - tau_ * lambda_values_loc[i][j];
    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i) {
      if (hyper_edge.node_descriptor[i] & (1<<6)) continue;  // static-only flag
      for (unsigned int j = 0; j < 2*space_dim; j++)
        if (hyper_edge.node_descriptor[i] & (1<<j))
          lambda_values_out[i][j] = 0.;
    }

    return lambda_values_out;
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
  std::array<lSol_float_t, 1U> errors(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    (void)lambda_values;

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    std::array<lSol_float_t,3> comps = {1,-1,-2};
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    lSol_float_t error = 0;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> r_old = hyper_edge.data.r_old;
      // loc_dof_to_glob_dof(hyper_edge.data.u_old, hyper_edge);

    for (unsigned int dim = 0; dim < 3; dim++) {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = u_old[i + dim * n_shape_fct_];
      error += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    for (unsigned int dim = 0; dim < 3; dim++) {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = r_old[i + dim * n_shape_fct_];
      error += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_phi, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    hy_check(std::isfinite(error), "error u_old" << u_old << " r_old " << r_old);

    return std::array<lSol_float_t, 1U>({error});
  }

  // Edge-local frame vector: idx == 0 → inner_normal(0),
  // idx == -k (k>=1) → outer_normal(k-1).
  // Convention: nonneg → inner, neg → outer.
  template <class hyEdgeT, typename float_t>
  Point<space_dim, float_t>
  edge_frame_vector(hyEdgeT& hyper_edge, int idx) const
  {
    if (idx >= 0)
      return (Point<space_dim, float_t>)hyper_edge.geometry.inner_normal(idx);
    else
      return (Point<space_dim, float_t>)hyper_edge.geometry.outer_normal(-idx - 1);
  }

  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    constexpr unsigned int n_pts    = Hypercube<hyEdge_dimT>::pow(sizeT);
    constexpr unsigned int n_fields = 6;  // n, m, u, r, v, s

   // Read stored fields directly (edge-local frame).
    // Order matches field index c: 0=n, 1=m, 2=u, 3=r, 4=v, 5=s.
    const SmallVec<space_dim*n_shape_fct_, lSol_float_t>* old_fields[n_fields] = {
      &hyper_edge.data.n_old,
      &hyper_edge.data.m_old,
      &hyper_edge.data.u_old,
      &hyper_edge.data.r_old,
      &hyper_edge.data.v_old,
      &hyper_edge.data.s_old,
    };

    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    SmallVec<static_cast<unsigned int>(sizeT), abscissa_float_t> helper(abscissas);

    // Edge-local-frame component values at each abscissa, per field.
    // shape: [field][local_component][pt]
    std::array<std::array<std::array<lSol_float_t, n_pts>, space_dim>, n_fields> point_vals{};

    for (unsigned int c = 0; c < n_fields; ++c)
      for (unsigned int dim = 0; dim < space_dim; ++dim) {
        for (unsigned int i = 0; i < coeffs.size(); ++i)
          coeffs[i] = (*old_fields[c])[dim * n_shape_fct_ + i];
        for (unsigned int pt = 0; pt < n_pts; ++pt)
          point_vals[c][dim][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
            coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT>>(pt, helper)
          );
      }

    std::array<std::array<lSol_float_t, n_pts>, system_dimension()> result{};

    static_assert(2 <= n_fields);

    // fields n m should be in local frame -> stay
    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int local = 0; local < space_dim; ++local)
        for (unsigned int q = 0; q < n_pts; ++q)
          result[c * space_dim + local][q] = point_vals[c][local][q];

    // fields u r should be in global frame -> need to transform
    for (unsigned int c = 2; c < n_fields; ++c)
      for (unsigned int local = 0; local < space_dim; ++local) {
        const int idx = -static_cast<int>(local);  // 0, -1, -2
        Point<space_dim, lSol_float_t> nv =
          edge_frame_vector<hyEdgeT, lSol_float_t>(hyper_edge, idx);
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          for (unsigned int q = 0; q < n_pts; ++q)
            result[c * space_dim + dim][q] += point_vals[c][local][q] * nv[dim];
      }

    return result;
  }

  template <class hyEdgeT>
  SmallVec<4*space_dim, lSol_float_t>
  get_extra_coeffs(hyEdgeT& hyper_edge) const
  {
    SmallVec<4 * space_dim, lSol_float_t> extra_coeffs(1.); // C_n, C_m, C_u, C_r

    if (!hyper_edge.geometry.has_extra_data()) return extra_coeffs;

    // 17 quantities associated with each fiber
    // column      header                  description
    //  0          mass,                   mass
    //  1, 2, 3    EA,kG_1A,kG_2A,         displacement stiffness
    //  4, 5, 6    G_xI_x,E_1I_1,E_2I_2,   rotation stiffness
    //  7, 8, 9    n_11,n_12,n_13,         normal 1
    // 10,11,12    n_21,n_22,n_23,         normal 2
    // 13,14       width1,width2,          widths in direction of normals
    // 15,16       fiber_id,fiber_edge_id  indicates which fiber this beam is part of
    //                                     -1 indicates no fiber, just virtual connection
    // material coefficients are given in the tangent,normal1,normal2 basis
    // so must be transformed into local basis chosen by HyperHDG, tangent coincides upto sign

    auto extra_data = hyper_edge.geometry.extra_data();
    hy_check(extra_data.size() == 17,
      "timowave expected 17 material coefficients, found " << extra_data.size());

    lSol_float_t mass = extra_data[0];
    SmallVec<space_dim, lSol_float_t> normal1 =
      std::array<lSol_float_t, space_dim>{{extra_data[7], extra_data[8], extra_data[9]}};
    SmallVec<space_dim, lSol_float_t> normal2 =
      std::array<lSol_float_t, space_dim>{{extra_data[10], extra_data[11], extra_data[12]}};
    SmallVec<space_dim, lSol_float_t> outer1 = hyper_edge.geometry.outer_normal(0);
    SmallVec<space_dim, lSol_float_t> outer2 = hyper_edge.geometry.outer_normal(1);

    auto length = hyper_edge.geometry.area();
    auto density = mass / (extra_data[13] * extra_data[14] * length);

    auto w1 = extra_data[13];
    auto w2 = extra_data[14];

    // approximation to second moment in direction of normal1, normal2
    // valid for rectangular cross section
    auto moment1 = w1*w1*w1 * w2 / 12;
    auto moment2 = w1 * w2*w2*w2 / 12;

    // displacement stiffness C_n
    extra_coeffs[0] = extra_data[1];
    extra_coeffs[1] = extra_data[2] * scalar_product(outer1, normal1) +
                      extra_data[3] * scalar_product(outer1, normal2);
    extra_coeffs[2] = extra_data[2] * scalar_product(outer2, normal1) +
                      extra_data[3] * scalar_product(outer2, normal2);

    // rotation stiffness C_m
    extra_coeffs[3] = extra_data[4];
    extra_coeffs[4] = extra_data[5] * scalar_product(outer1, normal1) +
                      extra_data[6] * scalar_product(outer1, normal2);
    extra_coeffs[5] = extra_data[5] * scalar_product(outer2, normal1) +
                      extra_data[6] * scalar_product(outer2, normal2);

    // displacement inertia C_u
    extra_coeffs[6] = mass / length;
    extra_coeffs[7] = mass / length;
    extra_coeffs[8] = mass / length;

    // rotation inertia
    extra_coeffs[9]  = density * (moment1 + moment2);
    extra_coeffs[10] = density * (moment1 * scalar_product(outer1, normal1) +
                                  moment2 * scalar_product(outer1, normal2));
    extra_coeffs[11] = density * (moment1 * scalar_product(outer2, normal1) +
                                  moment2 * scalar_product(outer2, normal2));

    for (unsigned int i = 0; i < extra_coeffs.size(); i++) {
      extra_coeffs[i] = std::abs(extra_coeffs[i]);
      hy_check(std::isfinite(extra_coeffs[i]), "get_extra finite coeffs " << extra_coeffs);
      hy_check(extra_coeffs[i] > 0, "get_extra zero coeffs " << extra_coeffs << " mass " << mass);
    }

    // hy_check(false, "extra_coeffs " << extra_coeffs);

    return extra_coeffs;
  }

  template <class hyEdgeT>
  void compute_fluxes(
    const std::array<std::array<lSol_float_t, 2*n_shape_bdr_*space_dim>, 2 * hyEdge_dimT>& lambda_values_loc,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& r_old = hyper_edge.data.r_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& n_old = hyper_edge.data.n_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& m_old = hyper_edge.data.m_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& v_old = hyper_edge.data.v_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& s_old = hyper_edge.data.s_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& flux_u = hyper_edge.data.flux_u;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& flux_r = hyper_edge.data.flux_r;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& flux_v = hyper_edge.data.flux_v;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& flux_s = hyper_edge.data.flux_s;

    SmallVec<4 * space_dim, lSol_float_t> extra_coeffs =
      get_extra_coeffs(hyper_edge);

    flux_u *= 0;
    flux_r *= 0;
    flux_v *= 0;
    flux_s *= 0;

    // NOTE: when we compute fluxes, -= for stuff from LHS, += for stuff from RHS
    //       finally is += to rhs

    for (unsigned int i = 0; i < n_shape_fct_; i++) {
      for (unsigned int j = 0; j < n_shape_fct_; j++) {
        SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec =
          integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
              decltype(hyEdgeT::geometry)>(j, i, hyper_edge.geometry);
        lSol_float_t bdr_int = 0;

        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          auto helper = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(i, j, face, hyper_edge.geometry);
          // grad_int_vec += helper * hyper_edge.geometry.local_normal(face);
          bdr_int += helper;
        }

        // std::cout << "---- compute fluxes" << std::endl;
        // std::cout << i << " " << j << "|" << bdr_int << " " << u_old[j] << std::endl;

        // NOTE: also need theta of old v with extra coeffs
        // NOTE: why no normal_int_vec here?
        for (unsigned int dim = 0; dim < space_dim; dim++) {
          flux_u[dim*n_shape_fct_ + i] += grad_int_vec[0]
            * n_old[dim*n_shape_fct_ +j] + tau_ * bdr_int * u_old[dim*n_shape_fct_+j];
          flux_r[dim*n_shape_fct_ + i] += grad_int_vec[0]
            * m_old[dim*n_shape_fct_ +j] + tau_ * bdr_int * r_old[dim*n_shape_fct_+j];
        }
      }

      for (unsigned int dim = 0; dim < space_dim; dim++) {
        flux_v[dim*n_shape_fct_ + i] += v_old[dim*n_shape_fct_+i] / extra_coeffs[2*space_dim+dim] * hyper_edge.geometry.area();
        flux_s[dim*n_shape_fct_ + i] += s_old[dim*n_shape_fct_+i] / extra_coeffs[3*space_dim+dim] * hyper_edge.geometry.area();
      }

      // Consider the cross product
      flux_r[2 * n_shape_fct_ + i] -= n_old[1 * n_shape_fct_ + i] * hyper_edge.geometry.area();
      flux_r[1 * n_shape_fct_ + i] += n_old[2 * n_shape_fct_ + i] * hyper_edge.geometry.area();

      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          for (unsigned int face = 0; face < 2 * hyEdge_dimT; face++) {
            flux_u[dim*n_shape_fct_+i] -= tau_*lambda_values_loc[face][j+dim*n_shape_bdr_]
              * integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(i,j, face, hyper_edge.geometry);
            flux_r[dim*n_shape_fct_+i] -= tau_*lambda_values_loc[face][j+(dim+space_dim)*n_shape_bdr_]
              * integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(i,j, face, hyper_edge.geometry);
          }
        }
      }
    }
  }

  template <class hyEdgeT>
  void set_data(
    const std::array<std::array<lSol_float_t, 2*n_shape_bdr_*space_dim>, 2 * hyEdge_dimT>& lambda_values_in,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    auto lambda_values = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    // std::cout << "  ---  set_data before" << std::endl;
    // std::cout << "v " << hyper_edge.data.v_old << std::endl;
    // std::cout << "s " << hyper_edge.data.s_old << std::endl;
    // std::cout << "lambda" << std::endl;
    // for (unsigned int i=0; i < lambda_values_in.size(); i++) {
    //   for (unsigned int j=0; j < lambda_values_in[i].size(); j++)
    //     std::cout << lambda_values[i][j] << " ";
    //   std::cout << std::endl;
    // }

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values, 1U, hyper_edge, time);

    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& r_old = hyper_edge.data.r_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& n_old = hyper_edge.data.n_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& m_old = hyper_edge.data.m_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& v_old = hyper_edge.data.v_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& s_old = hyper_edge.data.s_old;

    for (unsigned int i = 0; i < space_dim*n_shape_fct_; i++) {
      n_old[i] = coeffs[i+0*space_dim*n_shape_fct_];
      m_old[i] = coeffs[i+1*space_dim*n_shape_fct_];
      u_old[i] = coeffs[i+2*space_dim*n_shape_fct_];
      r_old[i] = coeffs[i+3*space_dim*n_shape_fct_];
      v_old[i] = coeffs[i+4*space_dim*n_shape_fct_];
      s_old[i] = coeffs[i+5*space_dim*n_shape_fct_];
    }

    compute_fluxes(lambda_values, hyper_edge, time);

    //std::cout << "----- set_data " << std::endl;
    //std::cout << "u " << hyper_edge.data.u_old << std::endl;
    //std::cout << "r " << hyper_edge.data.r_old << std::endl;
    //std::cout << "n " << hyper_edge.data.n_old << std::endl;
    //std::cout << "m " << hyper_edge.data.m_old << std::endl;
    //std::cout << "v " << hyper_edge.data.v_old << std::endl;
    //std::cout << "s " << hyper_edge.data.s_old << std::endl;
    //std::cout << "flux_u " << hyper_edge.data.flux_u << std::endl;
    //std::cout << "flux_r " << hyper_edge.data.flux_r << std::endl;
    //std::cout << "flux_v " << hyper_edge.data.flux_v << std::endl;
    //std::cout << "flux_s " << hyper_edge.data.flux_s << std::endl;
    //std::cout << "lambda_in" << std::endl;
    //for (unsigned int i=0; i < lambda_values_in.size(); i++) {
    //  for (unsigned int j=0; j < lambda_values_in[i].size(); j++)
    //    std::cout << lambda_values_in[i][j] << " ";
    //  std::cout << std::endl;
    //}
    //std::cout << "lambda " << std::endl;
    //for (unsigned int i=0; i < lambda_values.size(); i++) {
    //  for (unsigned int j=0; j < lambda_values[i].size(); j++)
    //    std::cout << lambda_values[i][j] << " ";
    //  std::cout << std::endl;
    //}

  }

  template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial(SmallMatT& lambda_values,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& r_old = hyper_edge.data.r_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& n_old = hyper_edge.data.n_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& m_old = hyper_edge.data.m_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& v_old = hyper_edge.data.v_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& s_old = hyper_edge.data.s_old;

    SmallVec<4 * space_dim, lSol_float_t> extra_coeffs(1.);

    // TODO: set v,s to something!!

    // first u then r in skeletal variables
    SmallVec<space_dim, lSol_float_t> helper;

    static_assert(lambda_values[0].size() == 2*space_dim);

    // Set skeltal variable!
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
        for (unsigned int j = 0; j < n_shape_bdr_; ++j) {
          helper = integrator::template integrate_bdrUni_psivecfunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::initial_u,
            Point<hyEdge_dimT, lSol_float_t>>(j, i, hyper_edge.geometry, time);
          for (unsigned int dim = 0; dim < space_dim; dim++)
            lambda_values[i][j+dim] = helper[dim];
           helper = integrator::template integrate_bdrUni_psivecfunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::initial_r,
            Point<hyEdge_dimT, lSol_float_t>>(j, i, hyper_edge.geometry, time);
          for (unsigned int dim = 0; dim < space_dim; dim++)
            lambda_values[i][j+space_dim+dim] = helper[dim];
        }
    }

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      for (unsigned int j = 0; j < 2*space_dim; ++j)
        if (hyper_edge.node_descriptor[i] & (1<<j))
            lambda_values[i][j] = 0.;

    auto lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    // set u, r, n, m old
    // NOTE: computing in global dofs
    // Define primary as L^2 projection!
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_u, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        u_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_v, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        v_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_s, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        s_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_r, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        r_old[i+dim*n_shape_fct_] = res[dim];
    }

    // std::cout << "----- make_initial (normals)" << std::endl;
    // std::cout << (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);
    // std::cout << (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(0);
    // std::cout << (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(1);

    // std::cout << "----- make_initial (glob)" << std::endl;
    // std::cout << "u " << hyper_edge.data.u_old;
    // std::cout << "v " << hyper_edge.data.v_old;
    // for (unsigned int i=0; i < lambda_values.size(); i++) {
    //   std::cout << "lambda " << i << "| ";
    //   for (unsigned int j=0; j < lambda_values[i].size(); j++)
    //     std::cout << lambda_values[i][j] << " ";
    //   std::cout << std::endl;
    // }

    // transform global dofs to edge dofs
    u_old = glob_dof_to_loc_dof(u_old, hyper_edge);
    r_old = glob_dof_to_loc_dof(r_old, hyper_edge);
    //n_old = glob_dof_to_loc_dof(n_old, hyper_edge);
    //m_old = glob_dof_to_loc_dof(m_old, hyper_edge);
    v_old = glob_dof_to_loc_dof(v_old, hyper_edge);
    s_old = glob_dof_to_loc_dof(s_old, hyper_edge);

    SmallVec<n_shape_fct_*6*space_dim, lSol_float_t> coeffs_old;
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int dim = 0; dim < space_dim; dim++) {
        coeffs_old[(2*space_dim+dim)*n_shape_fct_+i] = u_old[i+dim*n_shape_fct_];
      }
    }
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int dim = 0; dim < space_dim; dim++) {
        coeffs_old[(3*space_dim+dim)*n_shape_fct_+i] = r_old[i+dim*n_shape_fct_];
      }
    }

    auto mat = assemble_loc_matrix_a(hyper_edge, time);
    auto coefs_nm = mat*coeffs_old
      - assemble_rhs_from_lambda(lambda_values_loc, hyper_edge);
      //  - assemble_rhs_from_global_rhs(hyper_edge, time);

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        if (hyper_edge.node_descriptor[face]) {
          // u
          auto integrals1 = integrate_bdr_phivecfunccomp_beam<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::dirichlet_value_u>(i, face, {1,-1,-2}, hyper_edge.geometry, time);

          for (unsigned int comp = 0; comp < 3; comp++) {
            coefs_nm[(0 * space_dim + comp) * n_shape_fct_ + i] +=
              hyper_edge.geometry.local_normal(face).operator[](0) * integrals1[comp];
          }

          // phi
          integrals1 = integrate_bdr_phivecfunccomp_beam<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::dirichlet_value_phi>(i, face, {1,-1,-2}, hyper_edge.geometry, time);

          for (unsigned int comp = 0; comp < 3; comp++) {
            coefs_nm[(1 * space_dim + comp) * n_shape_fct_ + i] +=
              hyper_edge.geometry.local_normal(face).operator[](0) * integrals1[comp];
          }
        }
      }
    }

    // HACK: missing extra coefs -> function
    // though we dont use this function anyways if extra_coeffs != 1
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int d = 0; d < space_dim; ++d) {
        coefs_nm[i+(0*space_dim+d)*n_shape_fct_] *= -1;
      }
    }

    SmallSquareMat<n_shape_fct_> mmat;
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int j = 0; j < n_shape_fct_; ++j) {
        auto vol_integral = integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);
        mmat(i, j) = vol_integral;
      }
    }

    SmallVec<n_shape_fct_> temp;
    for (unsigned int d = 0; d < space_dim; d++) {
      for (unsigned int i = 0; i < n_shape_fct_; i++)
        temp[i] = coefs_nm[(0*space_dim+d)*n_shape_fct_+i];
      temp = temp / mmat;
      for (unsigned int i = 0; i < n_shape_fct_; i++)
        n_old[d*n_shape_fct_+i] = temp[i];
    }

    for (unsigned int d = 0; d < space_dim; d++) {
      for (unsigned int i = 0; i < n_shape_fct_; i++)
        temp[i] = coefs_nm[(1*space_dim+d)*n_shape_fct_+i];
      temp = temp / mmat;
      for (unsigned int i = 0; i < n_shape_fct_; i++)
        m_old[d*n_shape_fct_+i] = temp[i];
    }

    compute_fluxes(lambda_values_loc, hyper_edge, time);

    //std::cout << "----- make_initial" << std::endl;
    //std::cout << "u " << hyper_edge.data.u_old << std::endl;
    //std::cout << "r " << hyper_edge.data.r_old << std::endl;
    //std::cout << "n " << hyper_edge.data.n_old << std::endl;
    //std::cout << "m " << hyper_edge.data.m_old << std::endl;
    //std::cout << "v " << hyper_edge.data.v_old << std::endl;
    //std::cout << "s " << hyper_edge.data.s_old << std::endl;
    //for (unsigned int i=0; i < lambda_values.size(); i++) {
    //  std::cout << "lambda " << i << "| ";
    //  for (unsigned int j=0; j < lambda_values[i].size(); j++)
    //    std::cout << lambda_values_loc[i][j] << " ";
    //  std::cout << std::endl;
    //}
    //std::cout << "flux_u " << hyper_edge.data.flux_u << std::endl;
    //std::cout << "flux_r " << hyper_edge.data.flux_r << std::endl;
    //std::cout << "flux_v " << hyper_edge.data.flux_v << std::endl;
    //std::cout << "flux_s " << hyper_edge.data.flux_s << std::endl;

    return lambda_values;
  }

  template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial_from_static(SmallMatT& lambda_values_in,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    // HACK
    // solves static problem
    // with rhs from lambdas and global rhs but where we consider all dirichlet

    auto lambda_values = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    auto mat = assemble_loc_matrix_a(hyper_edge, time);
    SmallVec<n_loc_dofs_, lSol_float_t> rhs, coeffs;

    lSol_float_t integral;
    int comps[] = {1, -1, -2};
    static_assert(space_dim <= 3);

    // from lambda
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
          {
            rhs[(2 * space_dim + dim) * n_shape_fct_ + i] +=
              tau_ * lambda_values[face][j + dim] * integral;
            rhs[dim * n_shape_fct_ + i] -=
              hyper_edge.geometry.local_normal(face).operator[](0)
              * lambda_values[face][j + dim] * integral;
          }
        }

    // from global
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      for (unsigned int c = 0; c < space_dim; c++) {
        rhs[(2 * space_dim + c)* n_shape_fct_ + i] +=
          integrator::template integrate_vol_phivecfunccomp<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>
          >(i, comps[c], hyper_edge.geometry, 0.);

        rhs[(3 * space_dim + c)* n_shape_fct_ + i] +=
          integrator::template integrate_vol_phivecfunccomp<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>
          >(i, comps[c], hyper_edge.geometry, 0.);
      }
      // HACK: this only works if dirichlet_value_* returns zero in non-constrained directions
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face) {
        for (unsigned int c = 0; c < space_dim; c++) {
          if (hyper_edge.node_descriptor[face]) {
            integral = integrator::template integrate_bdr_phivecfunccomp<
              Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
              decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
              Point<hyEdge_dimT, lSol_float_t>>(i, face, comps[c], hyper_edge.geometry, 0.);
            rhs[(0 * space_dim + c) * n_shape_fct_ + i] -=
              hyper_edge.geometry.local_normal(face).operator[](0) * integral;
            rhs[(2 * space_dim + c) * n_shape_fct_ + i] += tau_ * integral;

            integral = integrator::template integrate_bdr_phivecfunccomp<
              Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
              decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
              Point<hyEdge_dimT, lSol_float_t>>(i, face, comps[c], hyper_edge.geometry, 0.);
            rhs[(1 * space_dim + c) * n_shape_fct_ + i] -=
              hyper_edge.geometry.local_normal(face).operator[](0) * integral;
            rhs[(3 * space_dim + c) * n_shape_fct_ + i] += tau_ * integral;
          }
        }
      }
    }

    coeffs = rhs / mat;

    //for (unsigned int i = 0; i < coeffs.size(); i++)
    //  hy_check(std::isfinite(coeffs[i]), "error coef? " << coeffs[i] << " coeffs " << coeffs << " rhs " << rhs << " mat " << mat);

    for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++) {
      hyper_edge.data.n_old[i] = coeffs[0*space_dim*n_shape_fct_+i];
      hyper_edge.data.m_old[i] = coeffs[1*space_dim*n_shape_fct_+i];
      hyper_edge.data.u_old[i] = coeffs[2*space_dim*n_shape_fct_+i];
      hyper_edge.data.r_old[i] = coeffs[3*space_dim*n_shape_fct_+i];
    }

    // rest is zero initialized

    // at the static-only dirichlet nodes the trace lambda is zero
    // at these dirichlet nodes, the trace lambda should be non zero for the wave problem,
    // hence we need to project the static bulk solution to the static-only dirichlet
    // trace lambda

    // Project bulk u_old, r_old onto trace lambda on bit-6 faces,
    // only in components where bit-j (j=0..2*space_dim-1) is unset.
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (!(hyper_edge.node_descriptor[face] & (1u << 6))) continue;

      for (unsigned int k = 0; k < n_shape_bdr_; ++k)
      {
        for (unsigned int d = 0; d < space_dim; ++d)
        {
            lSol_float_t num = 0;
            for (unsigned int i = 0; i < n_shape_fct_; ++i)
              num += integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                       i, k, face, hyper_edge.geometry)
                     * hyper_edge.data.u_old[d * n_shape_fct_ + i];
            lambda_values[face][k + d] = num;

            num = 0;
            for (unsigned int i = 0; i < n_shape_fct_; ++i)
              num += integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                       i, k, face, hyper_edge.geometry)
                     * hyper_edge.data.r_old[d * n_shape_fct_ + i];
            lambda_values[face][k + space_dim + d] = num;
        }
      }
    }

    compute_fluxes(lambda_values, hyper_edge, time);

    return lambda_values_in; // returns the input without changes
  }
};  // end of class LengtheningBernoulliBendingWave

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallSquareMat<
  TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_loc_matrix_a(hyEdgeT& hyper_edge, const lSol_float_t time) const
{
  // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;
  SmallVec<4 * space_dim, lSol_float_t> extra_coeffs = get_extra_coeffs(hyper_edge); // C_n, C_m, C_u, C_r

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);

      grad_int_vec =
        integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                       decltype(hyEdgeT::geometry)>(
          i, j, hyper_edge.geometry);

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        face_integral += helper;
        normal_int_vec += helper * hyper_edge.geometry.local_normal(face);
      }

      //TODO: mult by theta

      for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          vol_integral / extra_coeffs[dim];
        local_mat(dim * n_shape_fct_ + i, (2 * space_dim + dim) * n_shape_fct_ + j) -=
          grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          (normal_int_vec[0] - grad_int_vec[0]);
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                  (2 * space_dim + dim) * n_shape_fct_ + j) += tau_ * face_integral;

        local_mat((4 * space_dim + dim) * n_shape_fct_ + i,
                  (4 * space_dim + dim) * n_shape_fct_ + j) += vol_integral / extra_coeffs[2*space_dim+dim];
      }


      // Consider the cross product
      local_mat(2 * n_shape_fct_ + i, (3 * space_dim + 1) * n_shape_fct_ + j) -= vol_integral;
      local_mat(1 * n_shape_fct_ + i, (3 * space_dim + 2) * n_shape_fct_ + j) += vol_integral;
      local_mat((3 * space_dim + 2) * n_shape_fct_ + i, 1 * n_shape_fct_ + j) -= vol_integral;
      local_mat((3 * space_dim + 1) * n_shape_fct_ + i, 2 * n_shape_fct_ + j) += vol_integral;
    }
  }

  return local_mat;
}  // end of Diffusion::assemble_loc_matrix


template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallSquareMat<
  TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_loc_matrix(hyEdgeT& hyper_edge, const lSol_float_t time) const
{
  // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;
  SmallVec<4 * space_dim, lSol_float_t> extra_coeffs = get_extra_coeffs(hyper_edge); // C_n, C_m, C_u, C_r

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);

      grad_int_vec =
        integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                       decltype(hyEdgeT::geometry)>(
          i, j, hyper_edge.geometry);

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        face_integral += helper;
        normal_int_vec += helper * hyper_edge.geometry.local_normal(face);
      }

      //TODO: mult by theta

      for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          vol_integral / extra_coeffs[dim];
        local_mat(dim * n_shape_fct_ + i, (2 * space_dim + dim) * n_shape_fct_ + j) -=
          grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          theta_ * (normal_int_vec[0] - grad_int_vec[0]);
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                  (2 * space_dim + dim) * n_shape_fct_ + j) += theta_ * tau_ * face_integral;

        local_mat((4 * space_dim + dim) * n_shape_fct_ + i,
                  (4 * space_dim + dim) * n_shape_fct_ + j) += theta_ * vol_integral / extra_coeffs[2*space_dim+dim];
        local_mat((4 * space_dim + dim) * n_shape_fct_ + i,
                  (2 * space_dim + dim) * n_shape_fct_ + j) -= vol_integral / delta_t_;
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                  (4 * space_dim + dim) * n_shape_fct_ + j) += vol_integral / delta_t_;
      }


      // Consider the cross product
      local_mat(2 * n_shape_fct_ + i, (3 * space_dim + 1) * n_shape_fct_ + j) -= vol_integral;
      local_mat(1 * n_shape_fct_ + i, (3 * space_dim + 2) * n_shape_fct_ + j) += vol_integral;
      local_mat((3 * space_dim + 2) * n_shape_fct_ + i, 1 * n_shape_fct_ + j) -= theta_ * vol_integral;
      local_mat((3 * space_dim + 1) * n_shape_fct_ + i, 2 * n_shape_fct_ + j) += theta_ * vol_integral;
    }
  }

  return local_mat;
}  // end of Diffusion::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<
  TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_rhs_from_lambda(const SmallMatT& lambda_values, hyEdgeT& hyper_edge) const
{
  static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
                "Lambda values should have same floating point arithmetics as local solver!");
  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == 2 * space_dim * n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");

  // for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
  //   for (unsigned int j = 0; j < 1 * space_dim; ++j)
  //     std::cout << lambda_values[face][j] << " ";
  // std::cout << std::endl;

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
        {
          right_hand_side[(2 * space_dim + dim) * n_shape_fct_ + i] +=
            theta_ * tau_ * lambda_values[face][j + dim] * integral;
          right_hand_side[dim * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * lambda_values[face][j + dim] *
            integral;
        }
      }

  //std::cout << "-- rhs_from_lambda" << std::endl;
  //std::cout << right_hand_side << std::endl;
  return right_hand_side;
}  // end of Diffusion::assemble_rhs_from_lambda

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------


template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallVec<
  TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
// NOTE: dim unncessary
  assemble_rhs_from_global_rhs(hyEdgeT& hyper_edge, const lSol_float_t time) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  // constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  std::array<lSol_float_t, 3> integrals;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    // NOTE: it should probably not be *(space+dim) but *spac + dim???

    // distributed loads
    // f
    integrals = integrate_vol_phivecfunccomp_beam_avg<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n
      >(i, {1, -1, -2}, hyper_edge.geometry, time);
    for (unsigned int comp = 0; comp < 3; comp++)
      right_hand_side[(2*space_dim+comp) * n_shape_fct_ + i] = integrals[comp];

    // g
    integrals = integrate_vol_phivecfunccomp_beam_avg<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m
      >(i, {1, -1, -2}, hyper_edge.geometry, time);
    for (unsigned int comp = 0; comp < 3; comp++)
      right_hand_side[(3*space_dim+comp) * n_shape_fct_ + i] = integrals[comp];

    // NOTE: sign??
    // NOTE: think it should be subtracted here
    // NOTE: flux_* should be constructed with the sign as on the LHS, then it will be subtracted here

    // NOTE: it should probably not be *(space+dim) but *spac + dim???


    // NOTE: this fixes???
    for (unsigned int dim = 0; dim < space_dim; dim++) {
      right_hand_side[(2*space_dim+dim) * n_shape_fct_+i] -= (1-theta_)*hyper_edge.data.flux_u[dim*n_shape_fct_+i];
      right_hand_side[(3*space_dim+dim) * n_shape_fct_+i] -= (1-theta_)*hyper_edge.data.flux_r[dim*n_shape_fct_+i];
      right_hand_side[(4*space_dim+dim) * n_shape_fct_+i] -= (1-theta_)*hyper_edge.data.flux_v[dim*n_shape_fct_+i];
      right_hand_side[(5*space_dim+dim) * n_shape_fct_+i] -= (1-theta_)*hyper_edge.data.flux_s[dim*n_shape_fct_+i];
    }

    // std::cout << "  -- n" << std::endl;
    // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
    //   std::cout << right_hand_side[i+2*space_dim*n_shape_fct_] << " ";
    // std::cout << std::endl;

    // std::cout << "------ global_rhs" << std::endl;
    // dirichlet values
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (hyper_edge.node_descriptor[face] & (1<<6)) continue;
      if (hyper_edge.node_descriptor[face]) {
        // u
        auto integrals1 = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u>(i, face, {1,-1,-2}, hyper_edge.geometry, time);
        auto integrals2 = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u>(i, face, {1,-1,-2}, hyper_edge.geometry, time-delta_t_);

        for (unsigned int comp = 0; comp < 3; comp++) {
          right_hand_side[(0 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals1[comp];
          right_hand_side[(2 * space_dim + comp) * n_shape_fct_ + i] += tau_ * (theta_*integrals1[comp]+(1-theta_)*integrals2[comp]);
          // std::cout << i << " " << comp << "|" << theta_*integrals1[comp] << std::endl;
        }

        // phi
        integrals1 = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi>(i, face, {1,-1,-2}, hyper_edge.geometry, time);
        integrals2 = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi>(i, face, {1,-1,-2}, hyper_edge.geometry, time-delta_t_);

        for (unsigned int comp = 0; comp < 3; comp++) {
          right_hand_side[(1 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals1[comp];
          right_hand_side[(3 * space_dim + comp) * n_shape_fct_ + i] += tau_ * (theta_*integrals1[comp]+(1-theta_)*integrals2[comp]);
        }
      }
    }
  }

  // std::cout << "  -- n" << std::endl;
  // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  //   std::cout << right_hand_side[i+2*space_dim*n_shape_fct_] << " ";
  // std::cout << std::endl;

  // time derivatives
  for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  {
    right_hand_side[2 * space_dim * n_shape_fct_ + i] += hyper_edge.data.v_old[i] * hyper_edge.geometry.area() / delta_t_;
    right_hand_side[3 * space_dim * n_shape_fct_ + i] += hyper_edge.data.s_old[i] * hyper_edge.geometry.area() / delta_t_;
    right_hand_side[4 * space_dim * n_shape_fct_ + i] -= hyper_edge.data.u_old[i] * hyper_edge.geometry.area() / delta_t_;
    right_hand_side[5 * space_dim * n_shape_fct_ + i] -= hyper_edge.data.r_old[i] * hyper_edge.geometry.area() / delta_t_;
  }

  // std::cout << "-- rhs_from_global_rhs" << std::endl;
  // std::cout << "  -- u" << std::endl;
  // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  //   std::cout << right_hand_side[i] << " ";
  // std::cout << std::endl;
  // std::cout << "  -- r" << std::endl;
  // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  //   std::cout << right_hand_side[i+space_dim*n_shape_fct_] << " ";
  // std::cout << std::endl;
  // std::cout << "  -- n" << std::endl;
  // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  //   std::cout << right_hand_side[i+2*space_dim*n_shape_fct_] << " ";
  // std::cout << std::endl;
  // std::cout << "  -- m" << std::endl;
  // for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++)
  //   std::cout << right_hand_side[i+3*space_dim*n_shape_fct_] << " ";
  // std::cout << std::endl;

  // std::cout << "  -- v" << std::endl;
  // std::cout << hyper_edge.data.v_old << std::endl;
  // std::cout << "  -- s" << std::endl;
  // std::cout << hyper_edge.data.s_old << std::endl;
  // std::cout << std::endl;

  // std::cout << "  -- flux_u" << std::endl;
  // std::cout << hyper_edge.data.flux_u << std::endl;
  // std::cout << "  -- flux_r" << std::endl;
  // std::cout << hyper_edge.data.flux_r << std::endl;
  // std::cout << std::endl;


  //std::cout << "  -- rhs global" << std::endl;
  //std::cout << right_hand_side << std::endl;

  return right_hand_side;
}  // end of Bilaplacian::assemble_rhs_from_global_rhs

}  // namespace LocalSolver
