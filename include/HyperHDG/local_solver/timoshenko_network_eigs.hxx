#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <tuple>

namespace LocalSolver
{

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TimoschenkoBeamParametersDefault
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
   * \brief   Inverse diffusionbeam_network_bilaplacian.hxx coefficient in PDE as analytic function.
   ************************************************************************************************/
  // static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
  //                                              const param_float_t = 0.)
  // {
  //   return 1.;
  //   // return 1. / M_PI / M_PI;
  // }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return M_PI * (M_PI - 1.) * cos(M_PI * point[0]) * normal[2] *
             (point[1] == 0. && point[2] == 0.) +
           M_PI * cos(M_PI * point[1]) * (M_PI * normal[1] + normal[0]) *
             (point[0] == 0. && point[2] == 0.);
    // return -M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return (M_PI * M_PI - M_PI + 1.) * sin(M_PI * point[0]) * normal[1] *
             (point[1] == 0. && point[2] == 0.) +
           (M_PI * M_PI + 1.) * sin(M_PI * point[1]) * normal[2] *
             (point[0] == 0. && point[2] == 0.);
    // return  (M_PI * M_PI + 1.) * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2]
    // == 0.); return  -M_PI * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2] ==
    // 0.); return M_PI * M_PI * sin(M_PI * point[0]) * normal[0]; return M_PI * M_PI * sin(M_PI *
    // point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t initial_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return analytic_result_u(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return analytic_result_u(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t initial_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return analytic_result_phi(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return analytic_result_phi(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return 1.;
    return cos(M_PI * point[0]) * normal[2] + cos(M_PI * point[1]) * normal[1];
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return 1.;
    return sin(M_PI * point[0]) * normal[1] + sin(M_PI * point[1]) * normal[2];
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
};  // end of struct DiffusionParametersDefault

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TimoschenkoBeamParametersClamped
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
   * \brief   Inverse diffusionbeam_network_bilaplacian.hxx coefficient in PDE as analytic function.
   ************************************************************************************************/
  // static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
  //                                              const param_float_t = 0.)
  // {
  //   return 1.;
  //   // return 1. / M_PI / M_PI;
  // }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return 0;
    // return -M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return 0.;
    // return  (M_PI * M_PI + 1.) * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2]
    // == 0.); return  -M_PI * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2] ==
    // 0.); return M_PI * M_PI * sin(M_PI * point[0]) * normal[0]; return M_PI * M_PI * sin(M_PI *
    // point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return analytic_result_u(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return analytic_result_phi(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    // return 0.;
    return 5e-4 * ((point[0] > 1e-3) * normal[0] + (point[1] > .5e-3) * normal[1] + normal[2]);
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return 0.;
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
};  // end of struct DiffusionParametersDefault

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
          template <unsigned int, typename> typename parametersT = TimoschenkoBeamParametersDefault,
          typename lSol_float_t = double>
class TimoshenkoBeamEigs
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  struct data_type
  {
  };
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
  static constexpr unsigned int system_dimension() { return space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return space_dim; }

  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }

  template <class hyEdgeT>
  static constexpr std::array<bool, 2 * hyEdge_dimT>& dirichlet_nodes(std::array<bool, 2 * hyEdge_dimT>& input, const hyEdgeT& hyper_edge)
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      input[i] = is_dirichlet<parameters>(hyper_edge.node_descriptor[i]);
    return input;
  }

    template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial(SmallMatT& lambda_values,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          lambda_values[i][j] = 0.;
      else
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          lambda_values[i][(0 * space_dim + 0) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_u,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, 1, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);

        lambda_values[i][(0 * space_dim + 1) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_u,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, -1, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);

        lambda_values[i][(0 * space_dim + 2) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_u,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, -2, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);

        lambda_values[i][(1 * space_dim + 0) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_phi,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, 1, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);

        lambda_values[i][(1 * space_dim + 1) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_phi,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, -1, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);

        lambda_values[i][(1 * space_dim + 2) * n_shape_bdr_ + j] = integrator::template integrate_bdr_psivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::initial_phi,
          Point<hyEdge_dimT, lSol_float_t>>(j, i, -2, hyper_edge.geometry, 0.) / hyper_edge.geometry.face_area(i);
        }
    }

    return lambda_values;
  }

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
  static constexpr unsigned int n_loc_dofs_ = 4 * space_dim * n_shape_fct_;
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
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;

  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT>>,
    lSol_float_t>
    integrator;

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
  TimoshenkoBeamEigs(const constructor_value_type& tau = 1.) : tau_(tau) {}

  template <typename hyEdgeT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix(
    hyEdgeT& hyper_edge,
    const lSol_float_t eig) const;

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_lambda(
    const SmallMatT& lambda_values,
    hyEdgeT& hyper_edge) const;

  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs(
    hyEdgeT& hyper_edge,
    const unsigned int dim) const;

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
                                                                 hyEdgeT& hyper_edge,
                                                                 const lSol_float_t eig) const
  {
    try
    {
      SmallVec<n_loc_dofs_, lSol_float_t> rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
      return (rhs / assemble_loc_matrix(hyper_edge, eig)).data();
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
                              const lSol_float_t eig) const
  {
    hy_assert(lambda_values_in.size() == lambda_values_out.size() &&
                lambda_values_in.size() == 2 * hyEdge_dimT,
              "Both matrices must be of same size which corresponds to the number of faces!");
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      hy_assert(
        lambda_values_in[i].size() == lambda_values_out[i].size() &&
          lambda_values_in[i].size() == n_glob_dofs_per_node(),
        "Both matrices must be of same size which corresponds to the number of dofs per face!");

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        lambda_values_loc[i].fill(0.);

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, hyper_edge, eig);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = result(i, j) - tau_ * lambda_values_loc[i][j];

    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        lambda_values_out[i].fill(0.);

    return lambda_values_out;
  }


  template <class hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& jacobian_of_trace_to_flux(const SmallMatInT& lambda_values_in,
                                          SmallMatOutT& lambda_values_out,
                                          const lSol_float_t eig,
                                          const SmallMatInT& lambda_vals,
                                          const lSol_float_t eig_val,
                                          hyEdgeT& hyper_edge) const
  {
    hy_assert(lambda_values_in.size() == lambda_values_out.size() &&
                lambda_values_in.size() == 2 * hyEdge_dimT,
              "Both matrices must be of same size which corresponds to the number of faces!");
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      hy_assert(
        lambda_values_in[i].size() == lambda_values_out[i].size() &&
          lambda_values_in[i].size() == n_glob_dofs_per_node(),
        "Both matrices must be of same size which corresponds to the number of dofs per face!");

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        lambda_values_out[i].fill(0.);

    lambda_values_out = trace_to_flux(lambda_values_in, lambda_values_out, hyper_edge, eig_val);

    auto coeffs = solve_local_problem(lambda_vals, hyper_edge, eig_val);
    for (unsigned int i = 0; i < 2 * space_dim * n_shape_fct_; ++i)
      coeffs[i] = 0.;
    for (unsigned int i = 0; i < 2 * space_dim * n_shape_fct_; ++i)
      coeffs[2 * space_dim * n_shape_fct_ + i] *= eig * hyper_edge.geometry.area();

    coeffs = (SmallVec<coeffs.size(), lSol_float_t>(coeffs) /
              assemble_loc_matrix(hyper_edge, eig_val))
               .data();

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_out[i][j] += result(i, j);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        hy_assert(lambda_values_out[i][j] == lambda_values_out[i][j], "NaN");

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
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT, "Matrix must have appropriate size!");
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      hy_assert(lambda_values[i].size() == n_glob_dofs_per_node(),
                "Matrix must have appropriate size!");

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    auto lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      solve_local_problem(lambda_values_loc, 1U, hyper_edge, time);
    lSol_float_t error = 0.;
    std::array<lSol_float_t, n_shape_fct_> coeffs;

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 0) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, 1,
                                                                       hyper_edge.geometry, time);

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 1) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, -1,
                                                                       hyper_edge.geometry, time);

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 2) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, -2,
                                                                       hyper_edge.geometry, time);

    // error += len_beam.errors(lambda_values, hyper_edge, time)[0];
    // error += twist_beam.errors(lambda_values, hyper_edge, time)[0];
    // error += ben_beam.errors(lambda_values, hyper_edge, time)[0];
    return std::array<lSol_float_t, 1U>({error});
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
    // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    input_array_t lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);
    SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      solve_local_problem(lambda_values_loc, hyper_edge, time);
    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    SmallVec<static_cast<unsigned int>(sizeT), abscissa_float_t> helper(abscissas);

    // for(unsigned int i = 0; i < 2; ++i)
    //   for(unsigned int j = 0; j < 3; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl;

    // std::cout << coefficients;

    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
      point_vals, result;

    for (unsigned int dim = 0; dim < space_dim; ++dim)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[(2 * space_dim + dim) * n_shape_fct_ + i];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(sizeT); ++pt)
        point_vals[dim][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT>>(pt, helper));
    }

    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] = point_vals[0][q] * normal_vector[dim];

    normal_vector = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(0);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] += point_vals[1][q] * normal_vector[dim];

    normal_vector = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(1);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] += point_vals[2][q] * normal_vector[dim];

    return result;
  }
};  // end of class LengtheningBernoulliBendingBeam

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallSquareMat<
  TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_loc_matrix(hyEdgeT& hyper_edge, const lSol_float_t eig) const
{
  // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;
  SmallVec<2 * space_dim, lSol_float_t> extra_coeffs(1.);

  if (hyper_edge.geometry.has_extra_data())
  {
    auto extra_data = hyper_edge.geometry.extra_data();
    // (EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2):   6 structural constants
    // (n_11,n_12,n_13) : normal 1
    // (n_21,n_22,n_23) : normal 2
    SmallVec<space_dim, lSol_float_t> normal1 =
      std::array<lSol_float_t, space_dim>{{extra_data[6], extra_data[7], extra_data[8]}};
    SmallVec<space_dim, lSol_float_t> normal2 =
      std::array<lSol_float_t, space_dim>{{extra_data[9], extra_data[10], extra_data[11]}};
    SmallVec<space_dim, lSol_float_t> outer1 = hyper_edge.geometry.outer_normal(0);
    SmallVec<space_dim, lSol_float_t> outer2 = hyper_edge.geometry.outer_normal(1);

    extra_coeffs[0] = extra_data[0];
    extra_coeffs[1] = extra_data[1] * scalar_product(outer1, normal1) +
                      extra_data[2] * scalar_product(outer1, normal2);
    extra_coeffs[2] = extra_data[1] * scalar_product(outer2, normal1) +
                      extra_data[2] * scalar_product(outer2, normal2);
    extra_coeffs[3] = extra_data[3];
    extra_coeffs[4] = extra_data[4] * scalar_product(outer1, normal1) +
                      extra_data[5] * scalar_product(outer1, normal2);
    extra_coeffs[5] = extra_data[4] * scalar_product(outer2, normal1) +
                      extra_data[5] * scalar_product(outer2, normal2);

    extra_coeffs[0] *= 1e4;
    extra_coeffs[1] *= 1e4;
    extra_coeffs[2] *= 1e4;
    extra_coeffs[3] *= 1e12;
    extra_coeffs[4] *= 1e12;
    extra_coeffs[5] *= 1e12;
  }

  for (unsigned int i = 0; i < extra_coeffs.size(); ++i)
    extra_coeffs[i] = std::abs(extra_coeffs[i]);
  // extra_coeffs[i] = 1.;
  // hy_assert(false, extra_coeffs);

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

      for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          vol_integral / extra_coeffs[dim];
        local_mat(dim * n_shape_fct_ + i, (2 * space_dim + dim) * n_shape_fct_ + j) -=
          grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          normal_int_vec[0] - grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                  (2 * space_dim + dim) * n_shape_fct_ + j) += tau_ * face_integral;
      }

      // Consider the cross product
      local_mat(2 * n_shape_fct_ + i, (3 * space_dim + 1) * n_shape_fct_ + j) -= vol_integral;
      local_mat(1 * n_shape_fct_ + i, (3 * space_dim + 2) * n_shape_fct_ + j) += vol_integral;
      local_mat((3 * space_dim + 2) * n_shape_fct_ + i, 1 * n_shape_fct_ + j) -= vol_integral;
      local_mat((3 * space_dim + 1) * n_shape_fct_ + i, 2 * n_shape_fct_ + j) += vol_integral;

      for (unsigned int dim = 2 * space_dim; dim < 4 * space_dim; ++dim)
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) -= eig * vol_integral;
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
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<
  TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
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
            tau_ * lambda_values[face][j + dim] * integral;
          right_hand_side[dim * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * lambda_values[face][j + dim] *
            integral;
        }
      }

  // std::cout << right_hand_side << std::endl;
  return right_hand_side;
}  // end of Diffusion::assemble_rhs_from_lambda

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallVec<
  TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeamEigs<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_rhs_from_global_rhs(hyEdgeT& hyper_edge, const unsigned int dim) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  // constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[2 * space_dim * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, 1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(2 * space_dim + 1) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, -1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(2 * space_dim + 2) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, -2, hyper_edge.geometry,
                                                                         0.);

    right_hand_side[3 * space_dim * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, 1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(3 * space_dim + 1) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, -1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(3 * space_dim + 2) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, -2, hyper_edge.geometry,
                                                                         0.);
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[face]))
      {
        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, 1, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 0) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 0) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -1, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 1) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 1) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -2, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 2) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 2) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, 1, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 0) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 0) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -1, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 1) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 1) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -2, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 2) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 2) * n_shape_fct_ + i] += tau_ * integral;
      }
    }
  }

  return right_hand_side;
}  // end of Bilaplacian::assemble_rhs_from_global_rhs

}  // namespace LocalSolver
