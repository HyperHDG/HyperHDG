#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/ShapeFun1D.hxx>
#include <HyperHDG/QuadratureTensorial.hxx>
#include <HyperHDG/Hypercube.hxx>
#include <HyperHDG/DenseLA.hxx>

/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
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
template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t = double
>
class Diffusion_TensorialUniform
{
  public:
    /*!*********************************************************************************************
     * \brief Dimension of hyper edge type that this object solves on.
     * 
     * \todo  Why is this not just called dimension?
     *        -> E.g. in elasticity there are two important dimensions, the one of the hyperedge and
     *        the one of the space. Thus, elasticity can return both dimensions, while this class
     *        only returns the relevant hyperedge dimension.
     * \todo  The original brief referred to the internal variable only. It should be the other way
     *        round: this function is the main access to this number.
     *        -> I agree, this was not on purpose and I have to check for this in other classes!
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Decide whether gemetrical information is needed for local solver.
     * 
     * \retval  use_geom      True if geometrical information is used by local solver.
     **********************************************************************************************/
    static constexpr bool use_geometry() { return false; }
    /*!*********************************************************************************************
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
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return 1U; }
    
    static constexpr unsigned int system_dimension() { return 1U; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of quadrature points per spatial dimension.
     **********************************************************************************************/
    static constexpr unsigned int n_quads_1D_  = compute_n_quad_points<Gaussian>(quad_deg);
    /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1);
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node();
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = (hyEdge_dimT+1) * n_shape_fct_;
    /*!*********************************************************************************************
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
     **********************************************************************************************/
    static SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix ( const lSol_float_t tau );
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    /*!*********************************************************************************************
     * \brief   Local matrix for the local solver.
     **********************************************************************************************/
    const SmallSquareMat<n_loc_dofs_, lSol_float_t> loc_mat_;
    
    const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
    
    /*!*********************************************************************************************
     * \brief  Assemble local right hand for the local solver.
     *
     * The right hand side needs the values of the global degrees of freedom. Thus, it needs to be
     * constructed individually for every hyperedge.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs
    (const std::array<std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT>& lambda_values) const;
    
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem
    (const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
    {
      try { return (assemble_rhs(lambda_values) / loc_mat_).data(); }
      catch (LAPACKexception& exc)
      {
        hy_assert( 0 == 1 ,
                   exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        throw exc;
      }
    }
    /*!*********************************************************************************************
     * \brief   Evaluate primal variable at boundary.
     *
     * Function to evaluate primal variable of the solution. This function is needed to calculate
     * the local numerical fluxes.
     *
     * \param   coeffs        Coefficients of the local solution.
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > primal_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \param   coeffs        Coefficients of the local solution.
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, (hyEdge_dimT+1) * n_shape_fct_>& coeffs ) const;
  public:
    /*!*********************************************************************************************
     * \brief   Class is constructed using a single double indicating the penalty parameter.
     **********************************************************************************************/
    typedef lSol_float_t constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Constructor for local solver.
     *
     * \param   tau           Penalty parameter of HDG scheme.
     **********************************************************************************************/
    Diffusion_TensorialUniform(const constructor_value_type& tau = 1.)
    : tau_(tau), loc_mat_(assemble_loc_matrix(tau))  { } 
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \param   lambda_values Local part of vector x.
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    std::array< std::array<lSol_float_t, n_shape_bdr_>, 2 * hyEdge_dimT > numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
    {
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda_values);
      
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > 
        bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::system_dimension()>
    bulk_values
    (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values) const;

}; // end of class Diffusion_TensorialUniform


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF Diffusion_TensorialUniform
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
SmallSquareMat
<Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_, lSol_float_t>
Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::
assemble_loc_matrix ( const lSol_float_t tau )
{ 
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  lSol_float_t integral;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dimT>(i, j);
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        // Integral_element - nabla phi_i \vec phi_j dx 
        // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
        integral = integrator.template integrate_vol_Dphiphi<hyEdge_dimT>(i, j, dim);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) -= integral;
        local_mat(dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) -= integral;
    
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 1);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) += integral;
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 0);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) -= integral;
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
      }
    }
  }
  
  return local_mat;
} // end of Diffusion_TensorialUniform::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline SmallVec
<Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_, lSol_float_t>
Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::assemble_rhs
(const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
{
  lSol_float_t integral;

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;

  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j] * integral;
    
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j] * integral;
      }
    }
  }
  
  return right_hand_side;
} // end of Diffusion_TensorialUniform::assemble_rhs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::primal_at_boundary
( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
{
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        bdr_values[2*dim+0][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of Diffusion_TensorialUniform::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::dual_at_boundary
( const std::array<lSol_float_t, (hyEdge_dimT+1) * n_shape_fct_>& coeffs ) const
{
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        bdr_values[2*dim+0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of Diffusion_TensorialUniform::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t >
std::array
<
  std::array
  <
    lSol_float_t,
    Hypercube<hyEdge_dimT>::pow(sizeT)
  > ,
  Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::system_dimension()
>
Diffusion_TensorialUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::bulk_values
(const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values) const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients = solve_local_problem(lambda_values);

  std::array<std::array<lSol_float_t,Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()> values;
  std::array<unsigned int, hyEdge_dimT> dec_i, dec_q;
  lSol_float_t fct_value;
 
  std::array<unsigned int, poly_deg+1> poly_indices;
  for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
  std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
    values1D = shape_fct_eval<lSol_float_t,Legendre>(poly_indices, abscissas);
      
  for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  { 
    dec_i = integrator.template index_decompose<hyEdge_dimT>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = integrator.template index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of Diffusion_TensorialUniform::bulk_values





template < unsigned int space_dimT, typename param_float_t = double >
struct DiffusionParameters
{
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT,param_float_t>& pt)
  { return 1.;  }
  static param_float_t right_hand_side( const Point<space_dimT,param_float_t>& pt )
  { return 0.; }
};




/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
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
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters = DiffusionParameters<space_dimT,double>,
  typename lSol_float_t = double
>
class Diffusion
{
  public:
    /*!*********************************************************************************************
     * \brief Dimension of hyper edge type that this object solves on.
     * 
     * \todo  Why is this not just called dimension?
     *        -> E.g. in elasticity there are two important dimensions, the one of the hyperedge and
     *        the one of the space. Thus, elasticity can return both dimensions, while this class
     *        only returns the relevant hyperedge dimension.
     * \todo  The original brief referred to the internal variable only. It should be the other way
     *        round: this function is the main access to this number.
     *        -> I agree, this was not on purpose and I have to check for this in other classes!
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Decide whether gemetrical information is needed for local solver.
     * 
     * \retval  use_geom      True if geometrical information is used by local solver.
     **********************************************************************************************/
    static constexpr bool use_geometry() { return true; }
    /*!*********************************************************************************************
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
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return 1U; }
    
    static constexpr unsigned int system_dimension() { return 1U; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of quadrature points per spatial dimension.
     **********************************************************************************************/
    static constexpr unsigned int n_quads_1D_  = compute_n_quad_points<Gaussian>(quad_deg);
    /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1);
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node();
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = (hyEdge_dimT+1) * n_shape_fct_;
    /*!*********************************************************************************************
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
     **********************************************************************************************/
    template < typename GeomT > 
    inline SmallSquareMat<n_loc_dofs_, lSol_float_t>
    assemble_loc_matrix ( const lSol_float_t tau, GeomT& geom ) const;
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    
    const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
    
    /*!*********************************************************************************************
     * \brief  Assemble local right hand for the local solver.
     *
     * The right hand side needs the values of the global degrees of freedom. Thus, it needs to be
     * constructed individually for every hyperedge.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename GeomT >
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs
    ( const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values,
      GeomT& geom ) const;

    template < typename GeomT >
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs_from_global_rhs ( GeomT& geom ) const;
    
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename GeomT >
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem
    ( const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values,
      GeomT& geom ) const
    {
      try
      {
        return (assemble_rhs(lambda_values, geom) / assemble_loc_matrix(tau_, geom)).data();
      }
      catch (LAPACKexception& exc)
      {
        hy_assert( 0 == 1 ,
                   exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        throw exc;
      }
    }

    template < typename GeomT >
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem_rhs( GeomT& geom ) const
    {
      try
      {
        return (assemble_rhs_from_global_rhs(geom) / assemble_loc_matrix(tau_, geom)).data();
      }
      catch (LAPACKexception& exc)
      {
        hy_assert( 0 == 1 ,
                   exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        throw exc;
      }
    }
    /*!*********************************************************************************************
     * \brief   Evaluate primal variable at boundary.
     *
     * Function to evaluate primal variable of the solution. This function is needed to calculate
     * the local numerical fluxes.
     *
     * \param   coeffs        Coefficients of the local solution.
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename GeomT >
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dimT > primal_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs, GeomT& geom ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \param   coeffs        Coefficients of the local solution.
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename GeomT >
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, (hyEdge_dimT+1) * n_shape_fct_>& coeffs, GeomT& geom ) const;
  public:
    /*!*********************************************************************************************
     * \brief   Class is constructed using a single double indicating the penalty parameter.
     **********************************************************************************************/
    typedef lSol_float_t constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Constructor for local solver.
     *
     * \param   tau           Penalty parameter of HDG scheme.
     **********************************************************************************************/
    Diffusion(const constructor_value_type& tau = 1.) : tau_(tau)  { } 
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \param   lambda_values Local part of vector x.
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template <class GeomT>
    std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_lambda
    ( const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values,
      GeomT& geom ) const
    {
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda_values, geom);
      
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,geom)), duals(dual_at_boundary(coeffs,geom));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] 
            = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j] * geom.face_area(i);
       
      return bdr_values;
    }

    template <class GeomT>
    std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_rhs
    ( GeomT& geom ) const
    {
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem_rhs(geom);
      return dual_at_boundary(coeffs,geom);
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class GeomT>
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::system_dimension()
    >
    bulk_values
    ( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
      GeomT& geom ) const;

}; // end of class Diffusion


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF Diffusion
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename GeomT >
inline SmallSquareMat
<Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::n_loc_dofs_, lSol_float_t>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::
assemble_loc_matrix ( const lSol_float_t tau, GeomT& geom ) const
{ 
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator.template integrate_vol_phiphifunc
                        < GeomT, parameters::inverse_diffusion_coeff > (i, j, geom);
      // Integral_element - nabla phi_i \vec phi_j dx 
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
      grad_int_vec = integrator.template integrate_vol_nablaphiphi<GeomT>(i, j, geom);       

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator.template integrate_bdr_phiphi<GeomT>(i, j, face, geom);
        face_integral += helper;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          normal_int_vec[dim] += geom.hyEdge_dim_normal(face).operator[](dim) * helper; 
      }

      local_mat( hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) += tau * face_integral;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += vol_integral;
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) += normal_int_vec[dim];
      }
    }
  }
  
  return local_mat;
} // end of Diffusion::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename GeomT >
inline SmallVec
<Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::n_loc_dofs_, lSol_float_t>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::assemble_rhs
( const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values,
  GeomT& geom ) const
{
  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<GeomT>(i, j, face, geom);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * lambda_values[face][j] * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[dim * n_shape_fct_ + i]
            += geom.hyEdge_dim_normal(face).operator[](dim) * lambda_values[face][j] * integral; 
      }
  
  return right_hand_side;
} // end of Diffusion::assemble_rhs


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename GeomT >
inline SmallVec
<Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::n_loc_dofs_, lSol_float_t>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::
assemble_rhs_from_global_rhs ( GeomT& geom ) const
{
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    right_hand_side[hyEdge_dimT*n_shape_fct_ + i]
      = integrator.template integrate_vol_phifunc<GeomT,parameters::right_hand_side>(i, geom);
  
  return right_hand_side;
} // end of Diffusion::assemble_rhs_from_global_rhs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename GeomT >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::primal_at_boundary
( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs, GeomT& geom ) const
{
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        bdr_values[face][j] 
          += coeffs[hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<GeomT>(i, j, face, geom);
  
  return bdr_values;
} // end of Diffusion::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename GeomT >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::dual_at_boundary
( const std::array<lSol_float_t, (hyEdge_dimT+1) * n_shape_fct_>& coeffs, GeomT& geom ) const
{
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<GeomT>(i, j, face, geom);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          bdr_values[face][j] 
            += geom.hyEdge_dim_normal(face).operator[](dim) * integral
                 * coeffs[dim * n_shape_fct_ + i];
      }
  
  return bdr_values;
} // end of Diffusion::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, unsigned int poly_deg, unsigned int quad_deg,
  typename parameters, typename lSol_float_t
>
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, typename GeomT >
std::array
<
  std::array
  <
    lSol_float_t,
    Hypercube<hyEdge_dimT>::pow(sizeT)
  > ,
  Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::system_dimension()
>
Diffusion<hyEdge_dimT,space_dimT,poly_deg,quad_deg,parameters,lSol_float_t>::bulk_values
( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
  GeomT& geom ) const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients = solve_local_problem(lambda_values, geom);

  std::array<std::array<lSol_float_t,Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()> values;
  std::array<unsigned int, hyEdge_dimT> dec_i, dec_q;
  lSol_float_t fct_value;
 
  std::array<unsigned int, poly_deg+1> poly_indices;
  for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
  std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
    values1D = shape_fct_eval<lSol_float_t,Legendre>(poly_indices, abscissas);
      
  for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  { 
    dec_i = integrator.template index_decompose<hyEdge_dimT>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = integrator.template index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of Diffusion::bulk_values
