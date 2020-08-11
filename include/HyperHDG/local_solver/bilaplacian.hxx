#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/shape_fun_1d.hxx>
#include <HyperHDG/quadrature_tensorial.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/tensorial_shape_fun.hxx>
#include <algorithm>

/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo Update doxygen! Note that this may not depend on time (actually)!
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
class bilaplacian_uniform
{
  public:

    typedef struct empty_class {} data_type;

    typedef lSol_float_t solver_float_t;

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
    { return 2 * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return 1U; }
    
    static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = Hypercube<hyEdge_dimT>::pow(poly_deg + 1);
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1);
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = 2 * (hyEdge_dimT+1) * n_shape_fct_;
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
    (const std::array<std::array<lSol_float_t,2*n_shape_bdr_>, 2*hyEdge_dimT>& lambda_values) const;
    
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem
    (const std::array<std::array<lSol_float_t,2*n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
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
    inline std::array<std::array<lSol_float_t,2*n_shape_bdr_> , 2 * hyEdge_dimT> primal_at_boundary
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
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, 2 * (hyEdge_dimT+1) * n_shape_fct_>& coeffs ) const;
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
    bilaplacian_uniform(const constructor_value_type& tau = 1.)
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
    std::array< std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    ( const std::array< std::array<lSol_float_t, 2 * n_shape_bdr_>, 2 * hyEdge_dimT >&
      lambda_values, const lSol_float_t time = 0.
    )
    const
    {
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda_values);
      
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > 
        bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::system_dimension()>
    bulk_values
    ( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
      const lSol_float_t = 0.
    ) const;

  /*!*********************************************************************************************
   * \brief   Evaluate the function lambda on tensor product points on the boundary
   *
   * \tparam  absc_float_t  Floating type for the abscissa values.
   * \tparam  sizeT         Size of the array of array of abscissas.
   * \param   abscissas     Abscissas of the supporting points.
   * \param   lambda_values The values of the skeletal variable's coefficients.
   * \param   boundary_number number of the boundary on which to evaluate the function.
   **********************************************************************************************/
  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const;

}; // end of class bilaplacian_uniform


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF bilaplacian_uniform
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
SmallSquareMat
<bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_, lSol_float_t>
bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::
assemble_loc_matrix ( const lSol_float_t tau )
{ 
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  lSol_float_t integral;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dimT>(i, j);
      local_mat( hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j ) -= integral;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += integral;
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j ) += integral;
      }
      
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        // Integral_element - nabla phi_i \vec phi_j dx 
        // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
        integral = integrator.template integrate_vol_Dphiphi<hyEdge_dimT>(i, j, dim);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) -= integral;
        local_mat(dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) -= integral;
        local_mat(n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j)
          -= integral;
        local_mat(n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j)
          -= integral;
    
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 1);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) += integral;
        local_mat(n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j)
          += integral;  // Removing to enforce Neumann zero!
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
        local_mat(n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dimT>(i, j, 2 * dim + 0);
        local_mat(hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j) -= integral;
        local_mat(n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j)
           -= integral; // Removing to enforce Neumann zero!
        // and from the penalty in the lower right diagonal block
        local_mat(hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
        local_mat(n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j) 
          += tau * integral;
      }
    }
  }
  
  return local_mat;
} // end of BernoulliBendingBeam::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline SmallVec
<bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_, lSol_float_t>
bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::assemble_rhs
(const std::array< std::array<lSol_float_t, 2 * n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  lSol_float_t integral;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;

  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == 2 * n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        right_hand_side[dim*n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j] * integral;
        right_hand_side[n_dofs_lap + dim*n_shape_fct_ + i]
          += lambda_values[2*dim+0][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
          += tau_*lambda_values[2*dim+0][n_shape_bdr_ + j] * integral;
    
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j] * integral;
        right_hand_side[n_dofs_lap + dim*n_shape_fct_ + i]
          -= lambda_values[2*dim+1][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
          += tau_*lambda_values[2*dim+1][n_shape_bdr_ + j] * integral;
      }
    }
  }
  
  return right_hand_side;
} // end of BernoulliBendingBeam::assemble_rhs

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
    2*bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::primal_at_boundary
( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
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
        bdr_values[2*dim+0][n_shape_bdr_ + j] 
          += coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][j] += coeffs[hyEdge_dimT * n_shape_fct_ + i] * integral;
        bdr_values[2*dim+1][n_shape_bdr_ + j]
          += coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of bilaplacian_uniform::primal_at_boundary


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
    2*bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::dual_at_boundary
( const std::array<lSol_float_t, 2 * (hyEdge_dimT+1) * n_shape_fct_>& coeffs ) const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
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
        bdr_values[2*dim+0][n_shape_bdr_ + j]
          -= coeffs[n_dofs_lap + dim * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
        bdr_values[2*dim+1][n_shape_bdr_ + j]
          += coeffs[n_dofs_lap + dim * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of bilaplacian_uniform::dual_at_boundary


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
  bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::system_dimension()
>
bilaplacian_uniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>::bulk_values
(const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
 const lSol_float_t time ) const
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
    dec_i = index_decompose<hyEdge_dimT, poly_deg+1>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  return values;
}
template<unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t>
template<typename abscissa_float_t, std::size_t sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>, 1> bilaplacian_uniform<
    hyEdge_dimT,
    poly_deg,
    quad_deg,
    lSol_float_t>::lambda_values(const std::array<abscissa_float_t, sizeT> &abscissas,
                                   const input_array_t &lambda_values,
                                   const unsigned int boundary_number) const {
  return std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>, 1>();
}
// end of bilaplacian_uniform::bulk_values


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// GENERAL DIFFUSION PROBLEM AND RELATED CLASSES/STRUCTS
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct bilaplacian_parameters_default
{
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes {};
  static constexpr std::array<unsigned int, 0U> dirichlet_laplacian_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_laplacian_nodes {};
  
  static param_float_t inverse_bilaplacian_coefficient
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 1.; }
  static param_float_t right_hand_side
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t dirichlet_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t dirichlet_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t neumann_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t neumann_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t analytic_result
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
};
/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo    Update Doxygen! This may depend on time which has not been implemented so far!
 * 
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  - \nabla \cdot ( d \nabla u = f \quad \text{ in } \Omega, \qquad
 *  u = u_\textup D \quad \text{ on } \partial \Omega_\textup D}, \qquad
 *  - d \nabla u \cdot \nu = g_\textup N \quad \text{ on } \partial \Omega_\textup N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f$\f and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT    The dimension of the surrounding space.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains static parameter functions.
 *                        Defaults to above functions included in \c bilaplacianParametersDefault.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 * \tparam  space_dimTP   The dimension of the surrounding space.
 *                        Template parameter for the parameters which defaults to space_dimT.
 * \tparam  lSol_float_tP The floating point type calculations are executed in. Defaults to double.
 *                        Template parameter for the parameters which defaults to lSol_float_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT = bilaplacian_parameters_default,
  typename lSol_float_t = double
>
class bilaplacian
{
  public:
    
    typedef struct empty_class {} data_type;
    
    typedef lSol_float_t solver_float_t;

    // ---------------------------------------------------------------------------------------------
    // Public, static constexpr functions
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief Dimension of hyper edge type that this object solves on.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hypernode.
     **********************************************************************************************/
    static constexpr unsigned int node_value_dimension() { return 1U; }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }

  private:
  
    // ---------------------------------------------------------------------------------------------
    // Private, static constexpr functions
    // ---------------------------------------------------------------------------------------------

   /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1) / 2;
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node() / 2;
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = 2 * (hyEdge_dimT+1) * n_shape_fct_;

    template < typename parameters >
    static constexpr bool is_dirichlet( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                        node_type ) != parameters::dirichlet_nodes.end();
    }
    template < typename parameters >
    static constexpr bool is_dirichlet_laplacian( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_laplacian_nodes.begin(),
                        parameters::dirichlet_laplacian_nodes.end(),
                        node_type )
               != parameters::dirichlet_laplacian_nodes.end();
    }

    // ---------------------------------------------------------------------------------------------
    // Private, const members: Parameters and auxiliaries that help assembling matrices, etc.
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    /*!*********************************************************************************************
     * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
     **********************************************************************************************/
    const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
    
    // ---------------------------------------------------------------------------------------------
    // Private, internal functions for the local solver
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief  Assemble local matrix for the local solver.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   tau           Penalty parameter for HDG.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_mat       Matrix of the local solver.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline SmallSquareMat<n_loc_dofs_, lSol_float_t>
    assemble_loc_matrix
    ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t time ) const;
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from skeletal).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge 
    )  const;
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from global right-hand side).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs
    ( hyEdgeT& hyper_edge, const lSol_float_t time )  const;
    
    
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from global right-hand side).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t> assemble_rhs_from_coeffs
    ( const std::array< lSol_float_t, n_loc_dofs_ >&coeffs, hyEdgeT& hyper_edge )  const;
    
    
    /*!*********************************************************************************************
     * \brief  Solve local problem (with right-hand side from skeletal).
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      const unsigned int                                                            solution_type,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time
    )  const
    {
      try
      { 
        SmallVec<n_loc_dofs_, lSol_float_t> rhs;
        if (solution_type == 0)
          rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
        else if (solution_type == 1)
          rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge)
                  + assemble_rhs_from_global_rhs(hyper_edge, time);
        else  hy_assert( 0 == 1 , "This has not been implemented!" );
        return ( rhs / assemble_loc_matrix(tau_, hyper_edge, time) ).data();
      }
      catch (LAPACKexception& exc)
      {
        hy_assert( 0 == 1 ,
                   exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        throw exc;
      }
    }
    
    
    /*!*********************************************************************************************
     * \brief  Solve local problem (with right-hand side from skeletal).
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline std::array<lSol_float_t, n_loc_dofs_> solve_mass_problem
    ( 
      const std::array<lSol_float_t, n_loc_dofs_> & coeffs, hyEdgeT & hyper_edge,
      const lSol_float_t time
    )  const
    {
      try
      { 
        SmallVec<n_loc_dofs_, lSol_float_t> rhs = assemble_rhs_from_coeffs(coeffs, hyper_edge);
        return ( rhs / assemble_loc_matrix(tau_, hyper_edge, time) ).data();
      }
      catch (LAPACKexception& exc)
      {
        hy_assert( 0 == 1 ,
                   exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        throw exc;
      }
    }
    
    
        /*!*********************************************************************************************
     * \brief  Solve local problem (with right-hand side from skeletal).
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline std::array<lSol_float_t, n_loc_dofs_> solve_loc_prob_cor
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      const std::array<lSol_float_t, n_loc_dofs_> & coeffs,
      hyEdgeT                                                                   & hyper_edge,
      const lSol_float_t delta_time,
      const lSol_float_t time
    )  const
    {
      try
      { 
        SmallVec<n_loc_dofs_, lSol_float_t> rhs;
        rhs = delta_time * assemble_rhs_from_lambda(lambda_values, hyper_edge)
                + delta_time * assemble_rhs_from_global_rhs(hyper_edge,time)
                + assemble_rhs_from_coeffs(coeffs, hyper_edge);
        return ( rhs / (delta_time * assemble_loc_matrix(tau_, hyper_edge, time)) ).data();
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > primal_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs, hyEdgeT& hyper_edge ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_>& coeffs, hyEdgeT& hyper_edge ) const;
    
  public:
  
    // ---------------------------------------------------------------------------------------------
    // Public functions (and one typedef) to be utilized by external functions.
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief   Class is constructed using a single double indicating the penalty parameter.
     **********************************************************************************************/
    typedef lSol_float_t constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Constructor for local solver.
     *
     * \param   tau           Penalty parameter of HDG scheme.
     **********************************************************************************************/
    bilaplacian(const constructor_value_type& tau = 1.) : tau_(tau)  { } 
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 0U, hyper_edge, time);
      
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
              - tau_ * lambda_values[i][j < n_shape_bdr_ ? j+n_shape_bdr_ : j-n_shape_bdr_]
                     * hyper_edge.geometry.face_area(i);
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }
       
      return bdr_values;
    }

    /*!*********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * \todo    ALL!
   *
   * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
   * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
   * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
   * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
   *
   * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values Local part of vector x.
   * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
   * \retval  vecAx         Local part of vector A * x.
   **********************************************************************************************/
  template <class hyEdgeT>
  std::array<unsigned int, 2 * hyEdge_dimT> node_types(hyEdgeT& hyper_edge) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    
    std::array<unsigned int, 2 * hyEdge_dimT> result;
    result.fill(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))  result[i] = 1;

    return result;
  }
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_total
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 1U, hyper_edge, time);

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
                                 - tau_ * lambda_values[i][j] * hyper_edge.geometry.face_area(i);
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }
          
      return bdr_values;
    }
    
    
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_initial
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)  bdr_values[i][j] = 0.;
        else
        {
          for (unsigned int j = 0; j < lambda_values[i].size()/2; ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial >
                        (i, j, hyper_edge.geometry, time);
          for (unsigned int j = lambda_values[i].size()/2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial_laplace >
                        (i, j, hyper_edge.geometry, time);
        }
      }
          
      return bdr_values;
    }
    
    
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_mass
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 0U, hyper_edge, time);
      coeffs = solve_mass_problem(coeffs, hyper_edge, time);

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));

      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j];
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }

      return bdr_values;
    }
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > total_numerical_flux_mass
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 1U, hyper_edge, time);
      coeffs = solve_mass_problem(coeffs, hyper_edge, time);

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));

      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j];
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }

      return bdr_values;
    }
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    /*template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_mass
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 0U, hyper_edge, time);
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
      
      SmallVec<n_shape_fct_, lSol_float_t> u_coeffs, test_coeffs;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        u_coeffs[i] = coeffs[hyEdge_dimT*n_shape_fct_+i];

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > lambda_values_uni;
      for (unsigned int i = 0; i < lambda_values.size(); ++i)  lambda_values_uni[i].fill(0.);
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)  bdr_values[i][j] = 0.;
        else
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          {
            lambda_values_uni[i][j] = 1.;
            coeffs = solve_local_problem(lambda_values_uni, 0U, hyper_edge, time);
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
              test_coeffs[k] = coeffs[hyEdge_dimT*n_shape_fct_+k];
            bdr_values[i][j] = integrator.integrate_vol_phiphi
                                (u_coeffs.data(), test_coeffs.data(), hyper_edge.geometry);
            lambda_values_uni[i][j] = 0.;
          }

      return bdr_values;
    }*/
    
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    /*template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > total_numerical_flux_mass
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 1U, hyper_edge, time);
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
      
      SmallVec<n_shape_fct_, lSol_float_t> u_coeffs, test_coeffs;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        u_coeffs[i] = coeffs[hyEdge_dimT*n_shape_fct_+i];

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > lambda_values_uni;
      for (unsigned int i = 0; i < lambda_values.size(); ++i)  lambda_values_uni[i].fill(0.);
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)  bdr_values[i][j] = 0.;
        else
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          {
            lambda_values_uni[i][j] = 1.;
            coeffs = solve_local_problem(lambda_values_uni, 0U, hyper_edge, time);
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
              test_coeffs[k] = coeffs[hyEdge_dimT*n_shape_fct_+k];
            bdr_values[i][j] = integrator.integrate_vol_phiphi
                                (u_coeffs.data(), test_coeffs.data(), hyper_edge.geometry);
            lambda_values_uni[i][j] = 0.;
          }

      return bdr_values;
    }*/
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    lSol_float_t calc_L2_error_squared
    ( 
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                              & hy_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          hy_assert( lambda_values[i][j] == lambda_values[i][j] ,
                     "Lambda value wit index " << i << "," << j << " is NaN!" );
      
      std::array<lSol_float_t, n_loc_dofs_> coefficients
        = solve_local_problem(lambda_values, 1U, hy_edge, time);
      std::array< lSol_float_t, n_shape_fct_ > coeffs;
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + hyEdge_dimT * n_shape_fct_];
      
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        hy_assert( coeffs[i] == coeffs[i] , "The " << i << "-th coeff is NaN!" );
      
      lSol_float_t result = integrator.template integrate_vol_diffsquare_discana
          <decltype(hyEdgeT::geometry),parameters::analytic_result>(coeffs,hy_edge.geometry,time);
      hy_assert( result >= 0. , "The squared error must be non-negative, but was " << result );
      return result; 
    }
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    lSol_float_t calc_L2_error_squared_temp
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values_new,
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values_old,
      hyEdgeT                                                                   & hy_edge,
      const lSol_float_t delta_t,
      const lSol_float_t time
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      std::array<lSol_float_t, n_loc_dofs_> coeffs_new
        = solve_local_problem(lambda_values_new, 1U, hy_edge, time);      
      std::array<lSol_float_t, n_loc_dofs_> coeffs_old
        = solve_local_problem(lambda_values_old, 1U, hy_edge, time - delta_t);
      for (unsigned int i = 0; i < coeffs_old.size(); ++i)
        coeffs_old[i] = ( coeffs_old[i] - coeffs_new[i] );
      
      std::array<lSol_float_t, n_loc_dofs_> coefficients 
        = solve_loc_prob_cor(lambda_values_new, coeffs_old, hy_edge,  delta_t, time);
      std::array< lSol_float_t, n_shape_fct_ > coeffs;
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + hyEdge_dimT * n_shape_fct_];
      return integrator.template integrate_vol_diffsquare_discana
          <decltype(hyEdgeT::geometry),parameters::analytic_result>(coeffs,hy_edge.geometry,time);
    }

    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT >
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      bilaplacian
      < 
        hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t
      >::system_dimension()
    >
    bulk_values
    ( 
      const std::array<abscissa_float_t,sizeT>  & abscissas,
      const input_array_t                       & lambda_values,
      hyEdgeT                                   & hyper_edge,
      const lSol_float_t time = 0.
    )  const;
}; // end of class bilaplacian


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF bilaplacian
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallSquareMat
< bilaplacian <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_loc_matrix ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t time ) const
{ 
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, vol_func_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);
      vol_func_integral = integrator.template integrate_vol_phiphifunc
                          <decltype(hyEdgeT::geometry),parameters::inverse_bilaplacian_coefficient>
                          (i, j, hyper_edge.geometry, time);
      // Integral_element - nabla phi_i \vec phi_j dx 
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
      grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>
                       (i, j, hyper_edge.geometry);       

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>
                   (i, j, face, hyper_edge.geometry);
        face_integral += helper;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          normal_int_vec[dim] += hyper_edge.geometry.local_normal(face).operator[](dim) * helper;
      }

      local_mat( hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        -= vol_func_integral;

      local_mat( hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) += tau * face_integral;
      local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        += tau * face_integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += vol_integral;
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += vol_integral;
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          -= grad_int_vec[dim];
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
          -= grad_int_vec[dim];
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) += normal_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += normal_int_vec[dim];
      }
    }
  }
  
  return local_mat;
} // end of bilaplacian::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::assemble_rhs_from_lambda
( 
  const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
  hyEdgeT                                                                     & hyper_edge 
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == 2*n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * lambda_values[face][j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
          += tau_ * lambda_values[face][n_shape_bdr_ + j] * integral;
        
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          right_hand_side[dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][j] * integral;
          right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][n_shape_bdr_ + j] * integral;
        }
      }

  return right_hand_side;
} // end of bilaplacian::assemble_rhs_from_lambda


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_rhs_from_global_rhs ( hyEdgeT & hyper_edge, const lSol_float_t time )  const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
      = integrator.template integrate_vol_phifunc
          <decltype(hyEdgeT::geometry),parameters::right_hand_side>  (i, hyper_edge.geometry, time);
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[face]) )
      {
        integral = integrator.template integrate_bdr_phifunc
                     <decltype(hyEdgeT::geometry),parameters::dirichlet_value>
                     (i, face, hyper_edge.geometry, time);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[dim * n_shape_fct_ + i] 
            -= hyper_edge.geometry.local_normal(face).operator[](dim) * integral;
      }
      if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[face]) )
      {
        integral = integrator.template integrate_bdr_phifunc
                     <decltype(hyEdgeT::geometry),parameters::dirichlet_laplace_value>
                     (i, face, hyper_edge.geometry, time);
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i] += tau_ * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i] 
            -= hyper_edge.geometry.local_normal(face).operator[](dim) * integral;
      }
    }
  }
  return right_hand_side;
} // end of bilaplacian::assemble_rhs_from_global_rhs


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_coeffs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_rhs_from_coeffs
( 
  const std::array
  < 
    lSol_float_t,
    bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::n_loc_dofs_
  >& coeffs,
  hyEdgeT & hyper_edge
)  const
{
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      right_hand_side[ n_loc_dofs_/2 + hyEdge_dimT*n_shape_fct_ + i]
        += coeffs[ hyEdge_dimT*n_shape_fct_ + j] 
            * integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);
    }

  return right_hand_side;
} // end of Diffusion::assemble_rhs_from_coeffs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  < 
    lSol_float_t,
    2 * bilaplacian<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::primal_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_ >  & coeffs,
  hyEdgeT                                       & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        bdr_values[face][n_shape_bdr_ + j] 
          += coeffs[hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);
                  
        bdr_values[face][j] 
          += coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);          
      }
  
  return bdr_values;
} // end of bilaplacian::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    2 * bilaplacian<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::dual_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_>  & coeffs,
  hyEdgeT                                                         & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          bdr_values[face][n_shape_bdr_ + j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[dim * n_shape_fct_ + i];
          
          bdr_values[face][j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[n_dofs_lap + dim * n_shape_fct_ + i];
        }
      }
  
  return bdr_values;
} // end of bilaplacian::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, typename hyEdgeT >
std::array
<
  std::array < lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT) > ,
  bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::system_dimension()
>
bilaplacian < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::bulk_values
( 
  const std::array<abscissa_float_t,sizeT>  & abscissas,
  const input_array_t                       & lambda_values,
  hyEdgeT                                   & hyper_edge,
  const lSol_float_t time
)  const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients
    = solve_local_problem(lambda_values, 1U, hyper_edge, time);

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
    dec_i = index_decompose<hyEdge_dimT,poly_deg+1>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of bilaplacian::bulk_values















/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo    Update Doxygen! This may depend on time which has not been implemented so far!
 * 
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  - \nabla \cdot ( d \nabla u = f \quad \text{ in } \Omega, \qquad
 *  u = u_\textup D \quad \text{ on } \partial \Omega_\textup D}, \qquad
 *  - d \nabla u \cdot \nu = g_\textup N \quad \text{ on } \partial \Omega_\textup N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f$\f and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT    The dimension of the surrounding space.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains static parameter functions.
 *                        Defaults to above functions included in \c bilaplacian_parabParametersDefault.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 * \tparam  space_dimTP   The dimension of the surrounding space.
 *                        Template parameter for the parameters which defaults to space_dimT.
 * \tparam  lSol_float_tP The floating point type calculations are executed in. Defaults to double.
 *                        Template parameter for the parameters which defaults to lSol_float_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT,
  typename lSol_float_t = double
>
class bilaplacian_parab
{
  public:
    
    typedef lSol_float_t solver_float_t;

    // ---------------------------------------------------------------------------------------------
    // Public, static constexpr functions
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief Dimension of hyper edge type that this object solves on.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hypernode.
     **********************************************************************************************/
    static constexpr unsigned int node_value_dimension() { return 1U; }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }

  private:
  
    // ---------------------------------------------------------------------------------------------
    // Private, static constexpr functions
    // ---------------------------------------------------------------------------------------------

   /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1) / 2;
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node() / 2;
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = 2 * (hyEdge_dimT+1) * n_shape_fct_;

    template < typename parameters >
    static constexpr bool is_dirichlet( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                        node_type ) != parameters::dirichlet_nodes.end();
    }
    template < typename parameters >
    static constexpr bool is_dirichlet_laplacian( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_laplacian_nodes.begin(),
                        parameters::dirichlet_laplacian_nodes.end(),
                        node_type )
               != parameters::dirichlet_laplacian_nodes.end();
    }

    // ---------------------------------------------------------------------------------------------
    // Private, const members: Parameters and auxiliaries that help assembling matrices, etc.
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    
    const lSol_float_t theta_;

    const lSol_float_t delta_t_;
    
    /*!*********************************************************************************************
     * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
     **********************************************************************************************/
    const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
    
    // ---------------------------------------------------------------------------------------------
    // Private, internal functions for the local solver
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief  Assemble local matrix for the local solver.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   tau           Penalty parameter for HDG.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_mat       Matrix of the local solver.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline SmallSquareMat<n_loc_dofs_, lSol_float_t>
    assemble_loc_matrix
    ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t time ) const;
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from skeletal).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge 
    )  const;
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from global right-hand side).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs
    ( hyEdgeT& hyper_edge, const lSol_float_t time )  const;
    
    
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from global right-hand side).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t> assemble_rhs_from_coeffs
    ( const std::array< lSol_float_t, n_loc_dofs_ >&coeffs, hyEdgeT& hyper_edge )  const;
    
    
    /*!*********************************************************************************************
     * \brief  Solve local problem (with right-hand side from skeletal).
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      const unsigned int                                                            solution_type,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time
    )  const
    {
      try
      { 
        SmallVec<n_loc_dofs_, lSol_float_t> rhs;
        if (solution_type == 0)
          rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
        else if (solution_type == 1)
          rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge)
                  + assemble_rhs_from_global_rhs(hyper_edge, time);
        else  hy_assert( 0 == 1 , "This has not been implemented!" );
        return ( rhs / assemble_loc_matrix(tau_, hyper_edge, time) ).data();
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > primal_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs, hyEdgeT& hyper_edge ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_>& coeffs, hyEdgeT& hyper_edge ) const;
    
  public:
  
    // ---------------------------------------------------------------------------------------------
    // Public functions (and one typedef) to be utilized by external functions.
    // ---------------------------------------------------------------------------------------------
    
    typedef struct
    {
      std::array<lSol_float_t, n_loc_dofs_> coeffs;
    } data_type;
    
    /*!*********************************************************************************************
     * \brief   Class is constructed using a single double indicating the penalty parameter.
     **********************************************************************************************/
    typedef std::vector<lSol_float_t> constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Constructor for local solver.
     *
     * \param   tau           Penalty parameter of HDG scheme.
     **********************************************************************************************/
    bilaplacian_parab(const constructor_value_type& constru = std::vector(3,1.))
    : tau_(constru[0]), theta_(constru[1]), delta_t_(constru[2])  {}
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 0U, hyper_edge, time);
      
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
              - tau_ * lambda_values[i][j < n_shape_bdr_ ? j+n_shape_bdr_ : j-n_shape_bdr_]
                     * hyper_edge.geometry.face_area(i);
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }
       
      return bdr_values;
    }
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_total
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, 1U, hyper_edge, time);

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
                                 - tau_ * lambda_values[i][j] * hyper_edge.geometry.face_area(i);
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }
          
      return bdr_values;
    }
    
        
    /*!*********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
   * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
   * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
   * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
   *
   * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values Local part of vector x.
   * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
   * \retval  vecAx         Local part of vector A * x.
   **********************************************************************************************/
  template <class hyEdgeT>
  void set_data(
    const std::array<std::array<lSol_float_t,2*n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    std::array<lSol_float_t, n_loc_dofs_> coeffs =
      solve_local_problem(lambda_values, 1U, hyper_edge, time);

    hyper_edge.data.coeffs = coeffs;
  }
    
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_initial
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)  bdr_values[i][j] = 0.;
        else
        {
          for (unsigned int j = 0; j < lambda_values[i].size()/2; ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial >
                        (i, j, hyper_edge.geometry, time);
          for (unsigned int j = lambda_values[i].size()/2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial_laplace >
                        (i, j, hyper_edge.geometry, time);
        }
      }
      
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        hyper_edge.data.coeffs[hyEdge_dimT * n_shape_fct_ + i]
          = integrator.template integrate_volUni_phifunc<decltype(hyEdgeT::geometry), parameters::initial>(
                i, hyper_edge.geometry, time);
      
      
      return bdr_values;
    }
    
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    lSol_float_t calc_L2_error_squared
    ( 
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                              & hy_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          hy_assert( lambda_values[i][j] == lambda_values[i][j] ,
                     "Lambda value wit index " << i << "," << j << " is NaN!" );
      
      std::array<lSol_float_t, n_loc_dofs_> coefficients
        = solve_local_problem(lambda_values, 1U, hy_edge, time);
      std::array< lSol_float_t, n_shape_fct_ > coeffs;
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + hyEdge_dimT * n_shape_fct_];
      
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        hy_assert( coeffs[i] == coeffs[i] , "The " << i << "-th coeff is NaN!" );
      
      lSol_float_t result = integrator.template integrate_vol_diffsquare_discana
          <decltype(hyEdgeT::geometry),parameters::analytic_result>(coeffs,hy_edge.geometry,time);
      hy_assert( result >= 0. , "The squared error must be non-negative, but was " << result );
      return result; 
    }
    
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT >
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      bilaplacian_parab
      < 
        hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t
      >::system_dimension()
    >
    bulk_values
    ( 
      const std::array<abscissa_float_t,sizeT>  & abscissas,
      const input_array_t                       & lambda_values,
      hyEdgeT                                   & hyper_edge,
      const lSol_float_t time = 0.
    )  const;
}; // end of class bilaplacian_parab


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF bilaplacian_parab
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallSquareMat
< bilaplacian_parab <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_loc_matrix ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t time ) const
{ 
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, vol_func_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);
      vol_func_integral = integrator.template integrate_vol_phiphifunc
                          <decltype(hyEdgeT::geometry),parameters::inverse_bilaplacian_coefficient>
                          (i, j, hyper_edge.geometry, time);
      // Integral_element - nabla phi_i \vec phi_j dx 
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
      grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>
                       (i, j, hyper_edge.geometry);       

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>
                   (i, j, face, hyper_edge.geometry);
        face_integral += helper;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          normal_int_vec[dim] += hyper_edge.geometry.local_normal(face).operator[](dim) * helper;
      }

      local_mat( hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        -= vol_func_integral;

      local_mat( hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) += tau * face_integral;
      local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        += tau * theta_ * delta_t_ * face_integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += vol_integral;
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += vol_integral;
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          -= theta_ * delta_t_ * grad_int_vec[dim];
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
          -= grad_int_vec[dim];
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) += normal_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += theta_ * delta_t_ * normal_int_vec[dim];
      }
      
      local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) 
        += integrator.template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(i, j, hyper_edge.geometry);
      
    }
  }
  
  return local_mat;
} // end of bilaplacian_parab::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian_parab <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::assemble_rhs_from_lambda
( 
  const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
  hyEdgeT                                                                     & hyper_edge 
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == 2*n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * lambda_values[face][j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
          += tau_ * theta_ * delta_t_ * lambda_values[face][n_shape_bdr_ + j] * integral;
        
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          right_hand_side[dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][j] * integral;
          right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][n_shape_bdr_ + j] * integral;
        }
      }

  return right_hand_side;
} // end of bilaplacian_parab::assemble_rhs_from_lambda


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian_parab <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_rhs_from_global_rhs ( hyEdgeT & hyper_edge, const lSol_float_t time )  const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
      = theta_ * delta_t_ * integrator.template integrate_vol_phifunc
          <decltype(hyEdgeT::geometry),parameters::right_hand_side>  (i, hyper_edge.geometry, time)
        + (1 - theta_) * delta_t_ * integrator
          .template integrate_vol_phifunc<decltype(hyEdgeT::geometry), parameters::right_hand_side>(
            i, hyper_edge.geometry, time - delta_t_);
    
    // THIS HAS TO BE EXTENDED FOR THETA != 1
    
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[face]) )
      {
        integral = integrator.template integrate_bdr_phifunc
                     <decltype(hyEdgeT::geometry),parameters::dirichlet_value>
                     (i, face, hyper_edge.geometry, time);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[dim * n_shape_fct_ + i] 
            -= hyper_edge.geometry.local_normal(face).operator[](dim) * integral;
      }
      if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[face]) )
      {
        integral = integrator.template integrate_bdr_phifunc
                     <decltype(hyEdgeT::geometry),parameters::dirichlet_laplace_value>
                     (i, face, hyper_edge.geometry, time);
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i] += tau_ * theta_ * delta_t_ * integral;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i] 
            -= hyper_edge.geometry.local_normal(face).operator[](dim) * integral;
      }
    }
  }
  return right_hand_side + assemble_rhs_from_coeffs(hyper_edge.data.coeffs, hyper_edge);
} // end of bilaplacian_parab::assemble_rhs_from_global_rhs


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_coeffs
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian_parab <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_rhs_from_coeffs
( 
  const std::array
  < 
    lSol_float_t,
    bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::n_loc_dofs_
  >& coeffs,
  hyEdgeT & hyper_edge
)  const
{
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      right_hand_side[ n_loc_dofs_/2 + hyEdge_dimT*n_shape_fct_ + i]
        += coeffs[ hyEdge_dimT*n_shape_fct_ + j] 
            * integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);
    }

  return right_hand_side;
} // end of Diffusion::assemble_rhs_from_coeffs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  < 
    lSol_float_t,
    2 * bilaplacian_parab<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::primal_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_ >  & coeffs,
  hyEdgeT                                       & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        bdr_values[face][n_shape_bdr_ + j] 
          += coeffs[hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);
                  
        bdr_values[face][j] 
          += coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);          
      }
  
  return bdr_values;
} // end of bilaplacian_parab::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    2 * bilaplacian_parab<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::dual_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_>  & coeffs,
  hyEdgeT                                                         & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          bdr_values[face][n_shape_bdr_ + j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[dim * n_shape_fct_ + i];
          
          bdr_values[face][j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[n_dofs_lap + dim * n_shape_fct_ + i];
        }
      }
  
  return bdr_values;
} // end of bilaplacian_parab::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, typename hyEdgeT >
std::array
<
  std::array < lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT) > ,
  bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::system_dimension()
>
bilaplacian_parab < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::bulk_values
( 
  const std::array<abscissa_float_t,sizeT>  & abscissas,
  const input_array_t                       & lambda_values,
  hyEdgeT                                   & hyper_edge,
  const lSol_float_t time
)  const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients
    = solve_local_problem(lambda_values, 1U, hyper_edge, time);

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
    dec_i = index_decompose<hyEdge_dimT,poly_deg+1>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of bilaplacian_parab::bulk_values







/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo    Update Doxygen! This may depend on time which has not been implemented so far!
 * 
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  - \nabla \cdot ( d \nabla u = f \quad \text{ in } \Omega, \qquad
 *  u = u_\textup D \quad \text{ on } \partial \Omega_\textup D}, \qquad
 *  - d \nabla u \cdot \nu = g_\textup N \quad \text{ on } \partial \Omega_\textup N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f$\f and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT    The dimension of the surrounding space.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains static parameter functions.
 *                        Defaults to above functions included in \c bilaplacian_eigsParametersDefault.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 * \tparam  space_dimTP   The dimension of the surrounding space.
 *                        Template parameter for the parameters which defaults to space_dimT.
 * \tparam  lSol_float_tP The floating point type calculations are executed in. Defaults to double.
 *                        Template parameter for the parameters which defaults to lSol_float_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT,
  typename lSol_float_t = double
>
class bilaplacian_eigs
{
  public:
    
    typedef struct empty_class {} data_type;
    
    typedef lSol_float_t solver_float_t;

    // ---------------------------------------------------------------------------------------------
    // Public, static constexpr functions
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief Dimension of hyper edge type that this object solves on.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hypernode.
     **********************************************************************************************/
    static constexpr unsigned int node_value_dimension() { return 1U; }
    /*!*********************************************************************************************
     * \brief Dimension of of the solution evaluated with respect to a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }

  private:
  
    // ---------------------------------------------------------------------------------------------
    // Private, static constexpr functions
    // ---------------------------------------------------------------------------------------------

   /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1) / 2;
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node() / 2;
    /*!*********************************************************************************************
     * \brief   Number of (local) degrees of freedom per hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int n_loc_dofs_  = 2 * (hyEdge_dimT+1) * n_shape_fct_;

    template < typename parameters >
    static constexpr bool is_dirichlet( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                        node_type ) != parameters::dirichlet_nodes.end();
    }
    template < typename parameters >
    static constexpr bool is_dirichlet_laplacian( const unsigned int node_type )
    { 
      return std::find( parameters::dirichlet_laplacian_nodes.begin(),
                        parameters::dirichlet_laplacian_nodes.end(),
                        node_type )
               != parameters::dirichlet_laplacian_nodes.end();
    }

    // ---------------------------------------------------------------------------------------------
    // Private, const members: Parameters and auxiliaries that help assembling matrices, etc.
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    /*!*********************************************************************************************
     * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
     **********************************************************************************************/
    const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
    
    // ---------------------------------------------------------------------------------------------
    // Private, internal functions for the local solver
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief  Assemble local matrix for the local solver.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   tau           Penalty parameter for HDG.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_mat       Matrix of the local solver.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline SmallSquareMat<n_loc_dofs_, lSol_float_t>
    assemble_loc_matrix
    ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t eig ) const;
    /*!*********************************************************************************************
     * \brief  Assemble local right-hand for the local solver (from skeletal).
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_rhs       Local right hand side of the locasl solver.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge 
    )  const;  
     
    /*!*********************************************************************************************
     * \brief  Solve local problem (with right-hand side from skeletal).
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    template < typename hyEdgeT >  inline std::array<lSol_float_t, n_loc_dofs_> solve_local_problem
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t eig
    )  const
    {
      try
      { 
        SmallVec<n_loc_dofs_, lSol_float_t> rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
        return ( rhs / assemble_loc_matrix(tau_, hyper_edge, eig) ).data();
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
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > primal_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs, hyEdgeT& hyper_edge ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   coeffs        Coefficients of the local solution.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    template < typename hyEdgeT >
    inline std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT > dual_at_boundary
    ( const std::array<lSol_float_t, n_loc_dofs_>& coeffs, hyEdgeT& hyper_edge ) const;
    
  public:
  
    // ---------------------------------------------------------------------------------------------
    // Public functions (and one typedef) to be utilized by external functions.
    // ---------------------------------------------------------------------------------------------
    
    /*!*********************************************************************************************
     * \brief   Class is constructed using a single double indicating the penalty parameter.
     **********************************************************************************************/
    typedef lSol_float_t constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Constructor for local solver.
     *
     * \param   tau           Penalty parameter of HDG scheme.
     **********************************************************************************************/
    bilaplacian_eigs(const constructor_value_type& tau = 1.) : tau_(tau)  { } 
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_from_lambda
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t eig
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      std::array<lSol_float_t, n_loc_dofs_> coeffs
        = solve_local_problem(lambda_values, hyper_edge, eig);
      
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values,
        primals(primal_at_boundary(coeffs,hyper_edge)), duals(dual_at_boundary(coeffs,hyper_edge));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
              - tau_ * lambda_values[i][j < n_shape_bdr_ ? j+n_shape_bdr_ : j-n_shape_bdr_]
                     * hyper_edge.geometry.face_area(i);
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
            bdr_values[i][j] = 0.;
        if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = 0.;
      }
      
      return bdr_values;
    }
    
    /*!*********************************************************************************************
     * \todo    ALL
     **********************************************************************************************/    
    template <class hyEdgeT>
  std::array<unsigned int, 2 * hyEdge_dimT> node_types(hyEdgeT& hyper_edge) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    
    std::array<unsigned int, 2 * hyEdge_dimT> result;
    result.fill(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))  result[i] = 1;

    return result;
  }
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   lambda_values Local part of vector x.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vecAx         Local part of vector A * x.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > numerical_flux_initial
    ( 
      const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                                     & hyper_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

      std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
      {
        if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
          for (unsigned int j = 0; j < lambda_values[i].size(); ++j)  bdr_values[i][j] = 0.;
        else
        {
          for (unsigned int j = 0; j < lambda_values[i].size()/2; ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial >
                        (i, j, hyper_edge.geometry, time);
          for (unsigned int j = lambda_values[i].size()/2; j < lambda_values[i].size(); ++j)
            bdr_values[i][j] = integrator.template integrate_bdrUni_psifunc
              < decltype(hyEdgeT::geometry), parameters::initial_laplace >
                        (i, j, hyper_edge.geometry, time);
        }
      }
          
      return bdr_values;
    }
    
   /*!*********************************************************************************************
   * \brief   Evaluate local contribution to matrix--vector multiplication.
   *
   * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
   * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
   * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
   * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
   *
   * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values Local part of vector x.
   * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
   * \retval  vecAx         Local part of vector A * x.
   **********************************************************************************************/
  template <class hyEdgeT>
  std::array<std::array<lSol_float_t, 2*n_shape_bdr_>, 2 * hyEdge_dimT> numerical_flux_der(
    const std::array<std::array<lSol_float_t, 2*n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_values,
    const lSol_float_t eig,
    const std::array<std::array<lSol_float_t, 2*n_shape_bdr_>, 2 * hyEdge_dimT>& lambda_vals,
    const lSol_float_t eig_val,    
    hyEdgeT& hyper_edge) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    std::array<lSol_float_t, n_loc_dofs_> coeffs =
      solve_local_problem(lambda_values, hyper_edge, eig_val); // f_mu(eta,lambda_)

    std::array<std::array<lSol_float_t, 2*n_shape_bdr_>, 2 * hyEdge_dimT> bdr_values,
      primals(primal_at_boundary(coeffs, hyper_edge)), duals(dual_at_boundary(coeffs, hyper_edge));

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
      for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
        bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] 
            - tau_ * lambda_values[i][j < n_shape_bdr_ ? j+n_shape_bdr_ : j-n_shape_bdr_]
                   * hyper_edge.geometry.face_area(i);
      if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
        for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
          bdr_values[i][j] = 0.;
      if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
        for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = 0.;
    }
    
    coeffs = solve_local_problem(lambda_vals, hyper_edge, eig_val);
    
    // for (unsigned int i = 0; i < n_loc_dofs_ / 2 + hyEdge_dimT * n_shape_fct_; ++i)  coeffs[i] = 0.;
    // for (unsigned int i = 0; i < n_shape_fct_; ++i)  coeffs[n_loc_dofs_ / 2 + hyEdge_dimT * n_shape_fct_ + i] *= eig * hyper_edge.geometry.area();
    
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      coeffs[n_loc_dofs_ / 2 + hyEdge_dimT * n_shape_fct_ + i] = eig * coeffs[hyEdge_dimT * n_shape_fct_ + i] * hyper_edge.geometry.area();
    for (unsigned int i = 0; i < n_loc_dofs_ / 2 + hyEdge_dimT * n_shape_fct_; ++i)  coeffs[i] = 0.;
    
    coeffs = (SmallVec<coeffs.size(),lSol_float_t>(coeffs) / assemble_loc_matrix(tau_, hyper_edge, eig_val)).data();
    
    primals = primal_at_boundary(coeffs, hyper_edge);
    duals   = dual_at_boundary(coeffs, hyper_edge);
    
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
      for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
        bdr_values[i][j] += duals[i][j] + tau_ * primals[i][j];
      if ( is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) )
        for (unsigned int j = 0; j < lambda_values[i].size() / 2; ++j)
          bdr_values[i][j] = 0.;
      if ( is_dirichlet_laplacian<parameters>(hyper_edge.node_descriptor[i]) )
        for (unsigned int j = lambda_values[i].size() / 2; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = 0.;
    }
    
    return bdr_values;
  }
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    lSol_float_t calc_L2_error_squared
    ( 
      std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
      hyEdgeT                                                              & hy_edge,
      const lSol_float_t time = 0.
    )  const
    {
      using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
      
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          hy_assert( lambda_values[i][j] == lambda_values[i][j] ,
                     "Lambda value wit index " << i << "," << j << " is NaN!" );
      
      std::array<lSol_float_t, n_loc_dofs_> coefficients
        = solve_local_problem(lambda_values, hy_edge, time);
      std::array< lSol_float_t, n_shape_fct_ > coeffs;
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + hyEdge_dimT * n_shape_fct_];
      
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        hy_assert( coeffs[i] == coeffs[i] , "The " << i << "-th coeff is NaN!" );
      
      lSol_float_t result = integrator.template integrate_vol_diffsquare_discana
          <decltype(hyEdgeT::geometry),parameters::analytic_result>(coeffs,hy_edge.geometry,time);
      hy_assert( result >= 0. , "The squared error must be non-negative, but was " << result );
      return result; 
    }
    
    
    /*!*********************************************************************************************
     * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
     *
     * \tparam  absc_float_t  Floating type for the abscissa values.
     * \tparam  sizeT         Size of the array of array of abscissas.
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   abscissas     Abscissas of the supporting points.
     * \param   lambda_values The values of the skeletal variable's coefficients.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT >
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      bilaplacian_eigs
      < 
        hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t
      >::system_dimension()
    >
    bulk_values
    ( 
      const std::array<abscissa_float_t,sizeT>  & abscissas,
      const input_array_t                       & lambda_values,
      hyEdgeT                                   & hyper_edge,
      const lSol_float_t time = 0.
    )  const;
}; // end of class bilaplacian_eigs


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF bilaplacian_eigs
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallSquareMat
< bilaplacian_eigs <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::
assemble_loc_matrix ( const lSol_float_t tau, hyEdgeT& hyper_edge, const lSol_float_t eig ) const
{ 
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, vol_func_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator.template integrate_vol_phiphi(i, j, hyper_edge.geometry);
      vol_func_integral = integrator.template integrate_vol_phiphifunc
                          <decltype(hyEdgeT::geometry),parameters::inverse_bilaplacian_coefficient>
                          (i, j, hyper_edge.geometry, eig);
      // Integral_element - nabla phi_i \vec phi_j dx 
      // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
      grad_int_vec = integrator.template integrate_vol_nablaphiphi<decltype(hyEdgeT::geometry)>
                       (i, j, hyper_edge.geometry);       

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator.template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>
                   (i, j, face, hyper_edge.geometry);
        face_integral += helper;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          normal_int_vec[dim] += hyper_edge.geometry.local_normal(face).operator[](dim) * helper;
      }

      local_mat( hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        -= vol_func_integral;

      local_mat( hyEdge_dimT*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) += tau * face_integral;
      local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
        += tau * face_integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      { 
        local_mat( dim*n_shape_fct_+i , dim*n_shape_fct_+j ) += vol_integral;
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += vol_integral;
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( dim*n_shape_fct_+i , hyEdge_dimT*n_shape_fct_+j ) -= grad_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          -= grad_int_vec[dim];
        local_mat( n_dofs_lap + dim*n_shape_fct_+i , n_dofs_lap + hyEdge_dimT*n_shape_fct_+j )
          -= grad_int_vec[dim];
        
        local_mat( hyEdge_dimT*n_shape_fct_+i , dim*n_shape_fct_+j ) += normal_int_vec[dim];
        local_mat( n_dofs_lap + hyEdge_dimT*n_shape_fct_+i , n_dofs_lap + dim*n_shape_fct_+j )
          += normal_int_vec[dim];
      }
      
      local_mat(n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i, hyEdge_dimT * n_shape_fct_ + j) 
        -=  eig * integrator.template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(i, j, hyper_edge.geometry);
      
    }
  }
  
  return local_mat;
} // end of bilaplacian_eigs::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline SmallVec
< bilaplacian_eigs <hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_loc_dofs_, lSol_float_t >
bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::assemble_rhs_from_lambda
( 
  const std::array< std::array<lSol_float_t, 2*n_shape_bdr_>, 2*hyEdge_dimT > & lambda_values,
  hyEdgeT                                                                     & hyper_edge 
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  hy_assert( lambda_values.size() == 2 * hyEdge_dimT ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert( lambda_values[i].size() == 2*n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_ * lambda_values[face][j] * integral;
        right_hand_side[n_dofs_lap + hyEdge_dimT*n_shape_fct_ + i]
          += tau_ * lambda_values[face][n_shape_bdr_ + j] * integral;
        
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          right_hand_side[dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][j] * integral;
          right_hand_side[n_dofs_lap + dim * n_shape_fct_ + i]
            -= hyper_edge.geometry.local_normal(face).operator[](dim) 
                 * lambda_values[face][n_shape_bdr_ + j] * integral;
        }
      }

  return right_hand_side;
} // end of bilaplacian_eigs::assemble_rhs_from_lambda


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  < 
    lSol_float_t,
    2 * bilaplacian_eigs<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t>::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::primal_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_ >  & coeffs,
  hyEdgeT                                       & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        bdr_values[face][n_shape_bdr_ + j] 
          += coeffs[hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);
                  
        bdr_values[face][j] 
          += coeffs[n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i] 
              * integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                  (i, j, face, hyper_edge.geometry);          
      }
  
  return bdr_values;
} // end of bilaplacian_eigs::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename hyEdgeT >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    2 * bilaplacian_eigs<hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::n_shape_bdr_
  >,
  2 * hyEdge_dimT
>
bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::dual_at_boundary
( 
  const std::array<lSol_float_t, n_loc_dofs_>  & coeffs,
  hyEdgeT                                                         & hyper_edge
)  const
{
  constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator.template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>
                     (i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          bdr_values[face][n_shape_bdr_ + j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[dim * n_shape_fct_ + i];
          
          bdr_values[face][j] 
            += hyper_edge.geometry.local_normal(face).operator[](dim) * integral
                 * coeffs[n_dofs_lap + dim * n_shape_fct_ + i];
        }
      }
  
  return bdr_values;
} // end of bilaplacian_eigs::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int poly_deg, unsigned int quad_deg,
  template < unsigned int, typename >  typename parametersT, typename lSol_float_t
>
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, typename hyEdgeT >
std::array
<
  std::array < lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT) > ,
  bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::system_dimension()
>
bilaplacian_eigs < hyEdge_dimT,poly_deg,quad_deg,parametersT,lSol_float_t >::bulk_values
( 
  const std::array<abscissa_float_t,sizeT>  & abscissas,
  const input_array_t                       & lambda_values,
  hyEdgeT                                   & hyper_edge,
  const lSol_float_t time
)  const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients
    = solve_local_problem(lambda_values, hyper_edge, time);

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
    dec_i = index_decompose<hyEdge_dimT,poly_deg+1>(i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dimT>::pow(sizeT); ++q)
    {
      dec_q = index_decompose<hyEdge_dimT, abscissas.size()>(q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dimT; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of bilaplacian_eigs::bulk_values
