#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/quadrature/tensorial.hxx>
#include <HyperHDG/shape_function/shape_function.hxx>
#include <algorithm>
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
struct DiffusionParametersDefault
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes{};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
};  // end of struct DiffusionParametersDefault

/*!*************************************************************************************************
 * \brief   Local solver for parabolic diffusion equation on hypergraph.
 *
 * \note    Theta must be one, since only implicit Euler has been implemented at the moment.
 *
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  \partial_t u - \nabla \cdot ( d \nabla u ) = f \quad \text{ in } \Omega, \qquad
 *  u = u_\text D \quad \text{ on } \partial \Omega_\text D, \qquad
 *  - d \nabla u \cdot \nu = g_\text N \quad \text{ on } \partial \Omega_\text N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f\f$ and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains static parameter functions.
 *                        Defaults to above functions included in \c DiffusionParametersDefault.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = DiffusionParametersDefault,
          typename lSol_float_t = double>
class DiffusionParab
{
 public:
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
    return Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
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
   * \brief   Find out whether a node is of Dirichlet type.
   ************************************************************************************************/
  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }

  // -----------------------------------------------------------------------------------------------
  // Private, const members: Parameters and auxiliaries that help assembling matrices, etc.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;
  /*!***********************************************************************************************
   * \brief   Parameter theta that defines the one-step theta scheme.
   ************************************************************************************************/
  const lSol_float_t theta_;
  /*!***********************************************************************************************
   * \brief   Time step size.
   ************************************************************************************************/
  const lSol_float_t delta_t_;
  /*!***********************************************************************************************
   * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
   ************************************************************************************************/
  const Quadrature::Tensorial<
    Quadrature::Gaussian<quad_deg>,
    ShapeFunction<ShapeType::Tensorial<ShapeType::Legendre<poly_deg>, hyEdge_dimT> >,
    lSol_float_t>
    integrator;

  // -----------------------------------------------------------------------------------------------
  // Private, internal functions for the local solver
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Assemble local matrix for the local solver.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   tau           Penalty parameter.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Time, when the local matrix is evaluated.
   * \retval  loc_mat       Matrix of the local solver.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t>
  assemble_loc_matrix(const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t time) const;
  /*!***********************************************************************************************
   * \brief   Assemble local right-hand for the local solver (from skeletal).
   *
   * The right hand side needs the values of the global degrees of freedom. Note that we basically
   * follow the lines of
   *
   * B. Cockburn, J. Gopalakrishnan, and R. Lazarov.
   * Unified hybridization of discontinuous Galerkin, mixed, and continuous Galerkin methods for
   * second order elliptic problems. SIAM Journal on Numerical Analysis, 47(2):1319–1365, 2009
   *
   * and discriminate between local solvers with respect to the skeletal variable and with respect
   * to the global right-hand side. This assembles the local right-hand side with respect to the
   * skeletal variable.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of the \¢ lambda_values.
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \retval  loc_rhs       Local right hand side of the locasl solver.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_lambda(
    const SmallMatT& lambda_values,
    hyEdgeT& hyper_edge) const;
  /*!***********************************************************************************************
   * \brief   Assemble local right-hand for the local solver (from global right-hand side).
   *
   * Note that we basically follow the lines of
   *
   * B. Cockburn, J. Gopalakrishnan, and R. Lazarov.
   * Unified hybridization of discontinuous Galerkin, mixed, and continuous Galerkin methods for
   * second order elliptic problems. SIAM Journal on Numerical Analysis, 47(2):1319–1365, 2009
   *
   * and discriminate between local solvers with respect to the skeletal variable and with respect
   * to the global right-hand side. This assembles the local right-hand side with respect to the
   * global right-hand side. This function implicitly uses the parameters.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time, rhs is evaluated at
   * \retval  loc_rhs       Local right hand side of the locasl solver.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs(
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const;
  /*!***********************************************************************************************
   * \brief   Assemble local right-hand for the local solver (from function coefficients).
   *
   * Note that we basically follow the lines of
   *
   * B. Cockburn, J. Gopalakrishnan, and R. Lazarov.
   * Unified hybridization of discontinuous Galerkin, mixed, and continuous Galerkin methods for
   * second order elliptic problems. SIAM Journal on Numerical Analysis, 47(2):1319–1365, 2009
   *
   * and discriminate between local solvers with respect to the skeletal variable and with respect
   * to the global right-hand side. This assembles the local right-hand side with respect to the
   * global right-hand side. This function implicitly uses the parameters.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   coeffs        Local coefficients of bulk variable u.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \retval  loc_rhs       Local right hand side of the locasl solver.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_coeffs(
    const std::array<lSol_float_t, n_loc_dofs_>& coeffs,
    hyEdgeT& hyper_edge) const;
  /*!***********************************************************************************************
   * \brief   Solve local problem (with right-hand side from skeletal).
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of the \c lambda_values.
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \param   solution_type Type of local problem that is to be solved.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time the problem is solved.
   * \retval  loc_sol       Solution of the local problem.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem(const SmallMatT& lambda_values,
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
      return (rhs / assemble_loc_matrix(tau_, hyper_edge, time)).data();
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
   * \tparam  hyEdgeT      The geometry type / typename of the considered hyEdge's geometry.
   * \param   coeffs        Coefficients of the local solution.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> primal_at_boundary(
    const std::array<lSol_float_t, n_loc_dofs_>& coeffs,
    hyEdgeT& hyper_edge) const;
  /*!***********************************************************************************************
   * \brief   Evaluate dual variable at boundary.
   *
   * Function to evaluate dual variable of the solution. This function is needed to calculate the
   * local numerical fluxes.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   coeffs        Coefficients of the local solution.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> dual_at_boundary(
    const std::array<lSol_float_t, (hyEdge_dimT + 1) * n_shape_fct_>& coeffs,
    hyEdgeT& hyper_edge) const;

 public:
  // -----------------------------------------------------------------------------------------------
  // Public functions (and one typedef) to be utilized by external functions.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct
  {
    SmallVec<n_shape_fct_, lSol_float_t> u_old, flux_old;
  } data_type;
  /*!***********************************************************************************************
   *  \brief  Define type of node elements, especially with respect to nodal shape functions.
   ************************************************************************************************/
  typedef struct
  {
    typedef std::tuple<
      ShapeFunction<ShapeType::Tensorial<ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1> > >
      functions;
  } node_element;
  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef std::vector<lSol_float_t> constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   constru       Constructor object.
   ************************************************************************************************/
  DiffusionParab(const constructor_value_type& constru = std::vector(3, 1.))
  : tau_(constru[0]), theta_(constru[1]), delta_t_(constru[2])
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
   * \tparam  hyEdgeT             The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatInT         Data type of \c lambda_values_in.
   * \tparam  SmallMatOutT        Data type of \c lambda_values_out
   * \param   lambda_values_in    Local part of vector x.
   * \param   lambda_values_out   Local part that will be added to A * x.
   * \param   hyper_edge          The geometry of the considered hyperedge (of typename GeomT).
   * \param   time                Time at which time step ends.
   * \retval  vecAx               Local part of vector A * x.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_to_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time = 0.) const
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
    std::array<lSol_float_t, n_loc_dofs_> coeffs =
      solve_local_problem(lambda_values_in, 0U, hyper_edge, time);

    std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> primals(
      primal_at_boundary(coeffs, hyper_edge)),
      duals(dual_at_boundary(coeffs, hyper_edge));

    for (unsigned int i = 0; i < lambda_values_out.size(); ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] = 0.;
      else
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] =
            duals[i][j] + tau_ * primals[i][j] -
            tau_ * lambda_values_in[i][j] * hyper_edge.geometry.face_area(i);

    return lambda_values_out;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to residual.
   *
   * Execute residual evaluation y = A * x - b, where x represents the vector containing the
   * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
   * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
   * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
   *
   * \tparam  hyEdgeT             The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatInT         Data type of \c lambda_values_in.
   * \tparam  SmallMatOutT        Data type of \c lambda_values_out.
   * \param   lambda_values_in    Local part of vector x.
   * \param   lambda_values_out   Local part that will be added to A * x.
   * \param   hyper_edge          The geometry of the considered hyperedge (of typename GeomT).
   * \param   time                Time at which time step ends.
   * \retval  vecAx               Local part of vector A * x.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_and_data_to_flux(const SmallMatInT& lambda_values_in,
                                       SmallMatOutT& lambda_values_out,
                                       hyEdgeT& hyper_edge,
                                       const lSol_float_t time = 0.) const
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
    std::array<lSol_float_t, n_loc_dofs_> coeffs =
      solve_local_problem(lambda_values_in, 1U, hyper_edge, time);

    std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> primals(
      primal_at_boundary(coeffs, hyper_edge)),
      duals(dual_at_boundary(coeffs, hyper_edge));

    for (unsigned int i = 0; i < lambda_values_out.size(); ++i)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] = 0.;
      else
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] =
            duals[i][j] + tau_ * primals[i][j] -
            tau_ * lambda_values_in[i][j] * hyper_edge.geometry.face_area(i);
    }

    return lambda_values_out;
  }
  /*!***********************************************************************************************
   * \brief   Set data into data container to be used for next time step.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Time at which the old time step ended.
   ************************************************************************************************/
  template <class hyEdgeT>
  void set_data(
    const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    std::array<lSol_float_t, n_loc_dofs_> coeffs =
      solve_local_problem(lambda_values, 1U, hyper_edge, time);

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      hyper_edge.data.u_old[i] = coeffs[hyEdge_dimT * n_shape_fct_ + i];

    std::array<SmallVec<n_shape_fct_, lSol_float_t>, hyEdge_dimT> q_components;
    for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        q_components[dim][i] = coeffs[dim * n_shape_fct_ + i];

    // Fill flux_old!
    lSol_float_t helper;
    SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      hyper_edge.data.flux_old[i] =
        integrator
          .template integrate_vol_phifunc<decltype(hyEdgeT::geometry), parameters::right_hand_side>(
            i, hyper_edge.geometry, time);
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>(
          i, j, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          hyper_edge.data.flux_old[i] += q_components[dim][j] * grad_int_vec[dim];
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
            hyper_edge.data.flux_old[i] -= q_components[dim][j] *
                                           hyper_edge.geometry.local_normal(face).operator[](dim) *
                                           helper;
          hyper_edge.data.flux_old[i] -= tau_ * hyper_edge.data.u_old[j] * helper;
        }
      }
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
          hyper_edge.data.flux_old[i] +=
            tau_ * lambda_values[face][j] *
            integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              i, j, face, hyper_edge.geometry);
    }
  }
  /*!***********************************************************************************************
   * \brief   L2 project initial data to skeletal and fill data container.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data tyepe of \c lambda_values.
   * \param   lambda_values Local part of vector x.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Initial time of parabolic problem.
   * \retval  vecAx         L2 projected initial datum.
   ************************************************************************************************/
  template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial(SmallMatT& lambda_values,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    // Set skeltal variable!
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          lambda_values[i][j] = 0.;
      else
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          lambda_values[i][j] =
            integrator
              .template integrate_bdrUni_psifunc<decltype(hyEdgeT::geometry), parameters::initial>(
                j, i, hyper_edge.geometry, time);
    }

    // Define primary as L^2 projection!
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      hyper_edge.data.u_old[i] =
        integrator
          .template integrate_volUni_phifunc<decltype(hyEdgeT::geometry), parameters::initial>(
            i, hyper_edge.geometry, time);

    // Define dual as H^{1/2} projection!
    SmallSquareMat<n_shape_fct_, lSol_float_t> mass_flux;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        mass_flux(i, j) =
          integrator.template integrate_vol_phiphifunc<decltype(hyEdgeT::geometry),
                                                       parameters::inverse_diffusion_coeff>(
            i, j, hyper_edge.geometry, time);
    const SmallSquareMat<n_shape_fct_, lSol_float_t>& mass_mat = mass_flux;

    std::array<SmallVec<n_shape_fct_, lSol_float_t>, hyEdge_dimT> q_components;
    SmallVec<n_shape_fct_, lSol_float_t> local_rhs;

    for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
    {
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        local_rhs[i] =
          integrator
            .template integrate_vol_derphifunc<decltype(hyEdgeT::geometry), parameters::initial>(
              i, dim, hyper_edge.geometry, time);
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
          local_rhs[i] -=
            integrator
              .template integrate_bdr_phifunc<decltype(hyEdgeT::geometry), parameters::initial>(
                i, face, hyper_edge.geometry, time) *
            hyper_edge.geometry.local_normal(face).operator[](dim);
      }
      q_components[dim] = local_rhs / mass_mat;
    }

    // Fill flux_old!
    lSol_float_t helper;
    SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      hyper_edge.data.flux_old[i] =
        integrator
          .template integrate_vol_phifunc<decltype(hyEdgeT::geometry), parameters::right_hand_side>(
            i, hyper_edge.geometry, time);
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>(
          i, j, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          hyper_edge.data.flux_old[i] += q_components[dim][j] * grad_int_vec[dim];
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
            hyper_edge.data.flux_old[i] -= q_components[dim][j] *
                                           hyper_edge.geometry.local_normal(face).operator[](dim) *
                                           helper;
          hyper_edge.data.flux_old[i] -= tau_ * hyper_edge.data.u_old[j] * helper;
        }
      }
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
          hyper_edge.data.flux_old[i] +=
            tau_ * lambda_values[face][j] *
            integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              i, j, face, hyper_edge.geometry);
    }

    return lambda_values;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate squared local L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hy_edge           The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which error is evaluated.
   * \retval  err               Local squared L2 error.
   ************************************************************************************************/
  template <class hyEdgeT>
  lSol_float_t calc_L2_error_squared(const std::array<std::array<lSol_float_t, n_shape_bdr_>,
                                                      2 * hyEdge_dimT>& UNUSED(lambda_values),
                                     hyEdgeT& hy_edge,
                                     const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    return integrator.template integrate_vol_diffsquare_discana<decltype(hyEdgeT::geometry),
                                                                parameters::analytic_result>(
      hy_edge.data.u_old.data(), hy_edge.geometry, time);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  abscissa_float_t  Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Input array type.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which function is plotted.
   * \retval  func_vals         Array of function values.
   ************************************************************************************************/
  template <typename abscissa_float_t,
            std::size_t abscissas_sizeT,
            class input_array_t,
            class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>, system_dim>
  bulk_values(const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const;

};  // namespace LocalSolver

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF DiffusionParab
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallSquareMat<
  DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
  lSol_float_t>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::assemble_loc_matrix(
  const lSol_float_t tau,
  hyEdgeT& hyper_edge,
  const lSol_float_t time) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral =
        integrator.template integrate_vol_phiphifunc<decltype(hyEdgeT::geometry),
                                                     parameters::inverse_diffusion_coeff>(
          i, j, hyper_edge.geometry, time);
      // Integral_element - nabla phi_i \vec phi_j dx
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
      grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>(
        i, j, hyper_edge.geometry);

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        face_integral += helper;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          normal_int_vec[dim] += hyper_edge.geometry.local_normal(face).operator[](dim) * helper;
      }

      local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) +=
        tau * theta_ * delta_t_ * face_integral;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) += vol_integral;
        local_mat(dim * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) -= grad_int_vec[dim];
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) -=
          theta_ * delta_t_ * grad_int_vec[dim];
        local_mat(hyEdge_dimT * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          theta_ * delta_t_ * normal_int_vec[dim];
      }

      local_mat(hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) +=
        integrator.template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(i, j,
                                                                              hyper_edge.geometry);
    }
  }

  return local_mat;
}  // end of DiffusionParab::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<
  DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
  lSol_float_t>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_rhs_from_lambda(const SmallMatT& lambda_values, hyEdgeT& hyper_edge) const
{
  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
          tau_ * theta_ * delta_t_ * lambda_values[face][j] * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[dim * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](dim) * lambda_values[face][j] *
            integral;
      }

  return right_hand_side;
}  // end of DiffusionParab::assemble_rhs_from_lambda

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallVec<
  DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
  lSol_float_t>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_rhs_from_global_rhs(hyEdgeT& hyper_edge, const lSol_float_t time) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[hyEdge_dimT * n_shape_fct_ + i] =
      theta_ * delta_t_ *
      integrator
        .template integrate_vol_phifunc<decltype(hyEdgeT::geometry), parameters::right_hand_side>(
          i, hyper_edge.geometry, time);

    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (!is_dirichlet<parameters>(hyper_edge.node_descriptor[face]))
        continue;
      integral =
        integrator
          .template integrate_bdr_phifunc<decltype(hyEdgeT::geometry), parameters::dirichlet_value>(
            i, face, hyper_edge.geometry, time);
      right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
        tau_ * theta_ * delta_t_ * integral +
        tau_ * (1. - theta_) * delta_t_ *
          integrator.template integrate_bdr_phifunc<decltype(hyEdgeT::geometry),
                                                    parameters::dirichlet_value>(
            i, face, hyper_edge.geometry, time - delta_t_);
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        right_hand_side[dim * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](dim) * integral;
    }
  }

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
      hyper_edge.data.u_old[i] * hyper_edge.geometry.area() +
      delta_t_ * (1. - theta_) * hyper_edge.data.flux_old[i];

  return right_hand_side;
}  // end of DiffusionParab::assemble_rhs_from_global_rhs

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_coeffs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallVec<
  DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
  lSol_float_t>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::
  assemble_rhs_from_coeffs(
    const std::array<
      lSol_float_t,
      DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_>&
      coeffs,
    hyEdgeT& hyper_edge) const
{
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
      right_hand_side[hyEdge_dimT * n_shape_fct_ + i] +=
        coeffs[hyEdge_dimT * n_shape_fct_ + j] *
        integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);

  return right_hand_side;
}  // end of DiffusionParab::assemble_rhs_from_coeffs

// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline std::array<
  std::array<
    lSol_float_t,
    DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::primal_at_boundary(
  const std::array<lSol_float_t, n_loc_dofs_>& coeffs,
  hyEdgeT& hyper_edge) const
{
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)
    bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        bdr_values[face][j] +=
          coeffs[hyEdge_dimT * n_shape_fct_ + i] *
          integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);

  return bdr_values;
}  // end of DiffusionParab::primal_at_boundary

// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline std::array<
  std::array<
    lSol_float_t,
    DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_shape_bdr_>,
  2 * hyEdge_dimT>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::dual_at_boundary(
  const std::array<lSol_float_t, (hyEdge_dimT + 1) * n_shape_fct_>& coeffs,
  hyEdgeT& hyper_edge) const
{
  std::array<std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)
    bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          bdr_values[face][j] += hyper_edge.geometry.local_normal(face).operator[](dim) * integral *
                                 coeffs[dim * n_shape_fct_ + i];
      }

  return bdr_values;
}  // end of DiffusionParab::dual_at_boundary

// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename abscissa_float_t,
          std::size_t abscissas_sizeT,
          class input_array_t,
          typename hyEdgeT>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
           DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::bulk_values(
  const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
  const input_array_t& lambda_values,
  hyEdgeT& hyper_edge,
  const lSol_float_t time) const
{
  SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
    solve_local_problem(lambda_values, 1U, hyper_edge, time);
  SmallVec<n_shape_fct_, lSol_float_t> coeffs;
  SmallVec<static_cast<unsigned int>(abscissas_sizeT), abscissa_float_t> helper(abscissas);

  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
             DiffusionParab<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
    point_vals;

  for (unsigned int d = 0; d < system_dim; ++d)
  {
    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[d * n_shape_fct_ + i];
    for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
      point_vals[d][pt] = decltype(integrator)::shape_fun_t::template lin_comb_fct_val<float>(
        coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT> >(pt, helper));
  }

  return point_vals;
}
// end of DiffusionParab::bulk_values

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

}  // namespace LocalSolver
