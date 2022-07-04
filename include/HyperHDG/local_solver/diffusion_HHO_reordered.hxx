#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

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
 * \brief   Local solver diffusion equation on hypergraph.
 *
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  - \nabla \cdot ( d \nabla u ) = f \quad \text{ in } \Omega, \qquad
 *  u = u_\text D \quad \text{ on } \partial \Omega_\text D, \qquad
 *  - d \nabla u \cdot \nu = g_\text N \quad \text{ on } \partial \Omega_\text N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f\f$ and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * Moreover, this class provides the functions that approximate the condensed mass matrix. By this,
 * it is possible to approximate parabolic problems (in a very bad way) and to find good starting
 * values for the nonlinear eigenvalue problem.
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
class Diffusion_HHO
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
      TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1> > >
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
    return Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1); // q polynomials - maximum degree of each term is poly_deg, on faces (dimension 1 less than cell). thus (p+1)^(d-1)
  }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hyperedge.
   *
   *
   * \retval  n_dofs        Number of global degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_edge()
  {
    if (hyEdge_dimT==1)
      return Hypercube<hyEdge_dimT>::pow(poly_deg + 2);
    else
      return Hypercube<hyEdge_dimT>::pow(poly_deg + 2)-1; // polynomials with maximum degree polydeg+1 in each dimension on cells, thus (p+2)^d.  -1 (for 2D) because projection onto face should be nonzero, but for 3D?
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return  hyEdge_dimT + 1; }
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
  static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_edge(); 
  /*!***********************************************************************************************
   * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node(); 
  /*!***********************************************************************************************
   * \brief   Number of (local) degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = 2 * n_shape_fct_ - 1; // nshapefct for primal unknowns on the cell and nshapefct-1 for dual unknowns (since they are the gradient of the primal unknowns)
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
   * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
   ************************************************************************************************/
  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg+2>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg+1>, hyEdge_dimT> >,
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
   * \tparam  SmallMatT     The data type of the \c lambda values.
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
   * \brief   Solve local problem (with right-hand side from skeletal).
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of \c lambda_values.
   * \param   lambda_values Global degrees of freedom associated to the hyperedge.
   * \param   solution_type Type of local problem that is to be solved.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time the problem is solved.
   * \retval  loc_sol       Solution of the local problem.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> solve_local_problem(const SmallMatT& lambda_values,
                                                                 const unsigned int solution_type,
                                                                 hyEdgeT& hyper_edge,
                                                                 const lSol_float_t time) const
  {
    try
    {
      SmallVec<n_loc_dofs_, lSol_float_t> rhs;
      if (solution_type == 0)   // only homogeneous part (4.14)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
      else if (solution_type == 1)  // includes inhomogeneous part, i.e. (4.14) and (4.16)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge) +
              assemble_rhs_from_global_rhs(hyper_edge, time);
      else
        hy_assert(0 == 1, "This has not been implemented!");
      return rhs / assemble_loc_matrix(tau_, hyper_edge, time);
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
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \param   coeffs        Coefficients of the local solution.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> primal_at_boundary(
    const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
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
  inline SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> dual_at_boundary(
    const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
    hyEdgeT& hyper_edge) const;

 public:
  // -----------------------------------------------------------------------------------------------
  // Public functions (and one typedef) to be utilized by external functions.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef lSol_float_t constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  Diffusion_HHO(const constructor_value_type& tau = 1.) : tau_(tau) { }
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
   * \tparam  SmallMatOutT        Data type of \c lambda_values_out.
   *
   ********************************INPUT*******************************************************
   *
   * \param   lambda_values_in    Local part of vector x.
   * \param   lambda_values_out   Local part that will be added to A * x.
   * \param   hyper_edge          The geometry of the considered hyperedge (of typename GeomT).
   * \param   time                Time.
   *
   ********************************OUTPUT*******************************************************
   *
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
    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_in, 0U, hyper_edge, time);   // 4.14

    SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> primals(
      primal_at_boundary(coeffs, hyper_edge)),
      duals(dual_at_boundary(coeffs, hyper_edge));     

    for (unsigned int i = 0; i < lambda_values_out.size(); ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] = 0.;
      else
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
        {
          lambda_values_out[i][j] =
            duals(i, j) + tau_ * primals(i, j) -
            tau_ * lambda_values_in[i][j] * hyper_edge.geometry.face_area(i); // (the flux values are obtained from (4.14) -> which is why trace to flux, and then substituted to (4.18) before solving for \lambda)        
        }
        
    return lambda_values_out;
  }
  /*!***********************************************************************************************
   * \brief   Fill an array with 1 if the node is Dirichlet and 0 otherwise.
   ************************************************************************************************/
  template <typename hyEdgeT>
  std::array<unsigned int, 2 * hyEdge_dimT> node_types(hyEdgeT& hyper_edge) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    std::array<unsigned int, 2 * hyEdge_dimT> result;
    result.fill(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        result[i] = 1;

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local contribution to residual \f$ A x - b \f$.
   *
   * \tparam  hyEdgeT             The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatInT         Data type of \c lambda_values_in.
   * \tparam  SmallMatOutT        Data type of \c lambda_values_out.
   * \param   lambda_values_in    Local part of vector x.
   * \param   lambda_values_out   Local part that will be added to A * x.
   * \param   hyper_edge          The geometry of the considered hyperedge (of typename hyEdgeT).
   * \param   time                Time at which analytic functions are evaluated.
   * \retval  vecAx               Local part of vector A * x - b.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& residual_flux(const SmallMatInT& lambda_values_in,
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
    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_in, 1U, hyper_edge, time);  // 4.16

    SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> primals(
      primal_at_boundary(coeffs, hyper_edge)),
      duals(dual_at_boundary(coeffs, hyper_edge));
    
    for (unsigned int i = 0; i < lambda_values_out.size(); ++i)         // for each face
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)
          lambda_values_out[i][j] = 0.;
      else
        for (unsigned int j = 0; j < lambda_values_out[i].size(); ++j)   // for each basis function at the face
        {
          lambda_values_out[i][j] =
            duals(i, j) + tau_ * primals(i, j) -
            tau_ * lambda_values_in[i][j] * hyper_edge.geometry.face_area(i); // (the flux values are obtained from (4.14) -> which is why trace to flux, and then substituted to (4.18) before solving for \lambda)        //works because ONB?
        }
    }

    return lambda_values_out;
  }

  
  /*!***********************************************************************************************
   * \brief   Calculate squared local contribution of L2 error.
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     Data type of \c lambda_values.
   * \param   lambda_values The values of the skeletal variable's coefficients.
   * \param   hy_edge       The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Time at which anaytic functions are evaluated.
   * \retval  vec_b         Local part of vector b.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  std::array<lSol_float_t, 1U> errors(const SmallMatT& lambda_values,
                                      hyEdgeT& hy_edge,
                                      const lSol_float_t time = 0.) const
  {
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT, "Matrix must have appropriate size!");
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      hy_assert(lambda_values[i].size() == n_glob_dofs_per_node(),
                "Matrix must have appropriate size!");

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      solve_local_problem(lambda_values, 1U, hy_edge, time);
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i];

    return std::array<lSol_float_t, 1U>({integrator::template integrate_vol_diffsquare_discana<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result, Point<hyEdge_dimT, lSol_float_t> >(coeffs, hy_edge.geometry,
                                                                      time)});
  }

  /*!***********************************************************************************************
   * \brief   Evaluate local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  function_values   Function values at quadrature points.
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
};  // end of class Diffusion_HHO

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF Diffusion
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix (ok)
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallSquareMat<
  Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
  lSol_float_t>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::assemble_loc_matrix(
  const lSol_float_t tau,
  hyEdgeT& hyper_edge,
  const lSol_float_t time) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper, normal_ints_a, normal_ints_b, grad_int_vec, normal_int_vecA, normal_int_vecB;
  lSol_float_t beta_k, zeta_kphi_i;
//note: Here, first nshapefct unknowns correspond to that of U, the remaining for Q
  for (unsigned int i = 0; i < n_shape_fct_; ++i)   //ith test function
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j) // coeff of jth basis function
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      // note: the function here is 1/k, the inverse of the diffusion function
      vol_integral = integrator::template integrate_vol_nablaphinablaphifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::inverse_diffusion_coeff, Point<hyEdge_dimT, lSol_float_t> >(
        i, j, hyper_edge.geometry, time);
      if(i==0)
        hy_assert(vol_integral==0,"volume integral should be zero because gradient is zero, but value is" <<vol_integral);
      // Integral_element - nabla phi_i \vec phi_j dx
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
//      grad_int_vec =
//        integrator::template integrate_vol_nablaphinablaphifunc<
//        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
//        parameters::identity_coeff, Point<hyEdge_dimT, lSol_float_t> >(
//        i, j, hyper_edge.geometry, time);
      grad_int_vec =
        integrator::template integrate_vol_nablaphinablaphi<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        Point<hyEdge_dimT, lSol_float_t> >(
        i, j, hyper_edge.geometry, time);
      face_integral = 0.;
      normal_int_vecA = 0.;
      normal_int_vecB = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = 0.;
        for (unsigned int k = 0; k < n_shape_bdr_; ++k)
        {
          beta_k = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          j, k, face, hyper_edge.geometry) / hyper_edge.geometry.face_area(face); // coefficients of projection of phi_j onto the face
          zeta_kphi_i = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(i, k, face, hyper_edge.geometry); // integrate against test function
          helper += beta_k*zeta_kphi_i; // integral over the face
        }
//        if (i==2 && j==2 && face ==1)
//          hy_assert(1==0," The value of the integral at face" << face << " is " << helper );
        face_integral += helper;
//        normal_ints_a = integrator::template integrate_bdr_nablaphiphinufunc<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
//        parameters::identity_coeff, Point<hyEdge_dimT, lSol_float_t>>(i, j, face, hyper_edge.geometry, time);
        normal_ints_a = integrator::template integrate_bdr_nablaphiphinu<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        Point<hyEdge_dimT, lSol_float_t>>(i, j, face, hyper_edge.geometry);
        normal_int_vecA += normal_ints_a;
//        normal_ints_b = integrator::template integrate_bdr_nablaphiphinufunc<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
//        parameters::identity_coeff, Point<hyEdge_dimT, lSol_float_t>>(j, i, face, hyper_edge.geometry, time);
        normal_ints_b = integrator::template integrate_bdr_nablaphiphinu<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        Point<hyEdge_dimT, lSol_float_t>>(j, i, face, hyper_edge.geometry);
        normal_int_vecB += normal_ints_b;
        //if (j==1 && i==0 && face==1)
          //hy_assert(1==0, "gradient should be zero, but it is" <<normal_ints_b);
      }
//      if (j == n_shape_fct_-1 && i == n_shape_fct_-1 && hyEdge_dimT>1) 
//       local_mat(i,j) += tau;
//        hy_assert(1==0, "the value of the face integral here is" <<face_integral);
//      if (j==3 && i==0)
//        hy_assert(1==0, "normal ints should not be zero but the value for a and b are" <<normal_int_vecB << "and" <<normal_int_vecA);
// Order of equations: primal unknowns and test functions (4.14b) before dual (4.14a)
       // integral of tau* U_E v on the boundary in equation (4.14b)
       local_mat(i,j) += tau*face_integral;
        // equation (4.14b), integral of Q_E \cdot \nabla v - vQ \cdot n on the boundary
       local_mat(i,n_shape_fct_+j-1) -= (grad_int_vec - normal_int_vecB);
      // equation (4.14a), integral of Q_E \cdot p  
       local_mat(n_shape_fct_+i-1, n_shape_fct_+j-1) += vol_integral;      
     //equation (4.14a) integral of U_E div p
        local_mat(n_shape_fct_+i-1, j) += (grad_int_vec - normal_int_vecA);
      
    }
  }
//  if (hyEdge_dimT==2)
//    hy_assert(1==0, "this is the local matrix" << local_mat);
  return local_mat;
}  // end of Diffusion_HHO::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda (ok)
// -------------------------------------------------------------------------------------------------
// NOTE: This is independent of the RHS f
//
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
                lSol_float_t>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::assemble_rhs_from_lambda(
  const SmallMatT& lambda_values,
  hyEdgeT& hyper_edge) const
{
  static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
                "Lambda values should have same floating point arithmetics as local solver!");
  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral,integral_normal;
  lSol_float_t  time=0;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        integral_normal = integrator::template integrate_bdr_nablaphipsinufunc<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::identity_coeff, Point<hyEdge_dimT, lSol_float_t>>(
          i, j, face, hyper_edge.geometry, time);
          // rhs of equation (4.14b) note: lambda[f][j] is the coefficient of the jth basis function at face number f of the cell / hyperedge 
        right_hand_side[i] += tau_ * lambda_values[face][j] * integral;
        // rhs of equation (4.14a)
        right_hand_side[n_shape_fct_+i-1] -= lambda_values[face][j] * integral_normal;
      }
//  if(hyEdge_dimT==2)
//    hy_assert(0==1, "right hand side is " <<right_hand_side);
  return right_hand_side;
}  // end of Diffusion_HHO::assemble_rhs_from_lambda

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs (ok)
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t>
template <typename hyEdgeT>
inline SmallVec<Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_loc_dofs_,
                lSol_float_t>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::assemble_rhs_from_global_rhs(
  hyEdgeT& hyper_edge,
  const lSol_float_t time) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral, integral_normal;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[i] = integrator::template integrate_vol_phifunc<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::right_hand_side, Point<hyEdge_dimT, lSol_float_t> >(i, hyper_edge.geometry, time); // integral of fv in (4.16b)
    right_hand_side[n_shape_fct_+i-1] = 0.; 
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (!is_dirichlet<parameters>(hyper_edge.node_descriptor[face]))
        continue;
      integral = integrator::template integrate_bdr_phifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::dirichlet_value, Point<hyEdge_dimT, lSol_float_t> >(i, face,
                                                                        hyper_edge.geometry, time); // integral of u_Dv in (4.16b)
      integral_normal = integrator::template integrate_bdr_nablaphipsinufunc<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::dirichlet_value, Point<hyEdge_dimT, lSol_float_t>>(
          i, 0, face, hyper_edge.geometry, time); // integral of u_D p \cdot n in (4.16a)
//      hy_assert(integral_normal==0, "integral should be 0 since the dirichlet value is 0");
      right_hand_side[i] += tau_ * integral; 
      right_hand_side[n_shape_fct_+i-1] -= integral_normal; 
    }
  }

  return right_hand_side;
}  // end of Diffusion_HHO::assemble_rhs_from_global_rhs


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
inline SmallMat<2 * hyEdge_dimT,
                Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_shape_bdr_,
                lSol_float_t>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::primal_at_boundary(
  const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
  hyEdgeT& hyper_edge) const
{
  SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> bdr_values;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        bdr_values(face, j) +=
          coeffs[i] *
          integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);  // coefficient of the projection * basis function * mu_j

  return bdr_values;
}  // end of Diffusion_HHO::primal_at_boundary

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
inline SmallMat<2 * hyEdge_dimT,
                Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::n_shape_bdr_,
                lSol_float_t>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::dual_at_boundary(
  const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
  hyEdgeT& hyper_edge) const
{
using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallMat<2 * hyEdge_dimT, n_shape_bdr_, lSol_float_t> bdr_values;
  lSol_float_t integral;
  lSol_float_t time=0;

  for (unsigned int i = 1; i < n_shape_fct_; ++i) // starts with 1 because derivative of the basis function corresponding to i=0 is 0.
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator::template integrate_bdr_nablaphipsinufunc<Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::identity_coeff, Point<hyEdge_dimT, lSol_float_t>>(
          i, j, face, hyper_edge.geometry, time);
          bdr_values(face, j) += integral * coeffs[n_shape_fct_+i-1];   // integral over face of Q\cdot n * mu_j
      }

  return bdr_values;
}  // end of Diffusion_HHO::dual_at_boundary

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
           Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::bulk_values(
  const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
  const input_array_t& lambda_values,
  hyEdgeT& hyper_edge,
  const lSol_float_t time) const
{
  SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
    solve_local_problem(lambda_values, 1U, hyper_edge, time);
  SmallVec<n_shape_fct_, lSol_float_t> coeffs;
  Point<hyEdge_dimT> xyCoords;  
  
  SmallVec<static_cast<unsigned int>(abscissas_sizeT), abscissa_float_t> helper(abscissas);
  

  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
             Diffusion_HHO<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
    point_vals;

  for (unsigned int d = 0; d < system_dim; ++d)
  {
    for (unsigned int i = 0; i < coeffs.size(); ++i)
      	coeffs[i] = coefficients[i];
    for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
      if (d==system_dim-1)
        point_vals[d][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT> >(pt, helper));
      else
      {
        xyCoords = Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT>>(pt, helper);
        point_vals[d][pt] = xyCoords[d];
      }
  }

  return point_vals;
}
// end of Diffusion_HHO::bulk_values

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

}  // namespace LocalSolver
