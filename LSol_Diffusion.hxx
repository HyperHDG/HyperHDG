#ifndef LSOL_DIFFUSION_HXX
#define LSOL_DIFFUSION_HXX

#include <ShapeFun.hxx>
#include <Quadrature.hxx>
#include <FuncAndQuad.hxx>
#include <Hypercube.hxx>
#include <LapackWrapper.hxx>

/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * This class contains the local solver for Poisson's equation, i.e.,
 * \f[
 *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dim = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dim = 2) or volumes (\c hyEdge_dim = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t = double >
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
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
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
     *        problem (and therefore hyEdge_dim - dimensional). Thus, there is a discrimination
     *        between global and local amount per dofs in local solvers.
     * 
     *
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return Hypercube<hyEdge_dim-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return 1U; }
    
    static constexpr unsigned int system_dimension() { return 1U; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of quadrature points per spatial dimension.
     **********************************************************************************************/
    static constexpr unsigned int n_quads_1D_  = FuncQuad::compute_n_quad_points(quad_deg);
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
    static constexpr unsigned int n_loc_dofs_  = (hyEdge_dim+1) * n_shape_fct_;
    /*!*********************************************************************************************
     * \brief   Translate row and column indices to local index of entry in matrix.
     * 
     * Local \f$ n \times n \f$ matrices are encoded as arrays of size \f$n^2\f$. This function
     * translates a row and a column index into the index of the long array, where the corresponding
     * entry is located. Note that this is done column-wise (not row-wise as usually), to have the
     * correct format for LAPACK.
     *
     * The function is static inline, since it is used in the constructor's initializer list.
     *
     * \param   row           Row index of local mtatrix entry.
     * \param   column        Column index of local matrix entry.
     * \retval  index         Overall index of local matrix entry.
     **********************************************************************************************/
    static inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column)
    {
      hy_assert( 0 <= row && row < n_loc_dofs_ ,
                 "Row index must be >= 0 and smaller than total amount of local dofs." );
      hy_assert( 0 <= column && column < n_loc_dofs_ ,
                 "Column index must be >= 0 and smaller than total amount of local dofs." );

      return column * n_loc_dofs_ + row;  // Transposed for LAPACK
    }
    /*!*********************************************************************************************
     * \brief   Decompose index of local tensorial shape function into its (spatial) components.
     *
     * The function is static inline, since it is used in the constructor's initializer list.
     *
     * \param   index         Index of tensorial shape function.
     * \param   range         Amount of 1D shape functions (the tensorial is made up from).
     * \param   decomposition Array to be filled with the decomposition.
     * \retval  decomposition Array filled with the decomposition.
     **********************************************************************************************/
    template<unsigned int dimT> static inline void index_decompose ( unsigned int index, 
      unsigned int range, std::array<unsigned int, std::max(dimT,1U)>& decomposition )
    {
      if ( decomposition.size() == 1 )  decomposition[0] = index;
      else
      {
        for (unsigned int dim = 0; dim < dimT; ++dim)
        {
          decomposition[dim] = index % range;
          index /= range;
        }
      }
    }
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
    static std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_ >
    assemble_loc_matrix ( const lSol_float_t tau );
    /*!*********************************************************************************************
     * \brief   (Globally constant) penalty parameter for HDG scheme.
     **********************************************************************************************/
    const lSol_float_t tau_;
    /*!*********************************************************************************************
     * \brief   Quadrature weights per spatial dimension.
     **********************************************************************************************/
    const std::array<lSol_float_t, n_quads_1D_> q_weights_;
    /*!*********************************************************************************************
     * \brief   Trial functions evaluated at quadrature points (per spatial dimensions).
     **********************************************************************************************/
    const std::array< std::array<lSol_float_t, n_quads_1D_ > , poly_deg + 1 > trial_;
    /*!*********************************************************************************************
     * \brief   Trial functions evaluated at boundaries {0,1} (per spatial dimension).
     **********************************************************************************************/
    const std::array< std::array<lSol_float_t, 2 > , poly_deg + 1 > trial_bdr_;
    /*!*********************************************************************************************
     * \brief   Local matrix for the local solver.
     **********************************************************************************************/
    const std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_ > loc_mat_;
    
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
    inline std::array<lSol_float_t, n_loc_dofs_ > assemble_rhs
    (const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dim >& lambda_values) const;
    
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem
    (const std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
    {
      // A copy of loc_mat_ is created, since LAPACK will destroy the matrix values.
      std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_> local_matrix = loc_mat_;
      // The local right hand side is assembled (and will also be destroyed by LAPACK).
      std::array<lSol_float_t, n_loc_dofs_> right_hand_side = assemble_rhs(lambda_values);
      // LAPACK solves local_matix * return_value = right_hand_side.
      try { return lapack_solve<n_loc_dofs_>(local_matrix, right_hand_side); }
      catch (LASolveException& exc)
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
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > primal_at_boundary
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
    inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > dual_at_boundary
    ( const std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_>& coeffs ) const;
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
    Diffusion_TensorialUniform(const constructor_value_type& tau)
    : tau_(tau), q_weights_(FuncQuad::quad_weights<quad_deg,lSol_float_t>()),
      trial_(FuncQuad::shape_fcts_at_quad_points<poly_deg, quad_deg, lSol_float_t>()),
      trial_bdr_(FuncQuad::shape_fcts_at_bdrs<poly_deg, lSol_float_t>()),
      loc_mat_(assemble_loc_matrix(tau))
    { } 
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
    std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
    {
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda_values);
      
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > 
        bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(sizeT)>,
      Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()>
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
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
std::array
< 
  lSol_float_t,
  Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_ *
  Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_
>
Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::
assemble_loc_matrix ( const lSol_float_t tau )
{ 
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  lSol_float_t integral;
  
  std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_> local_mat;
  local_mat.fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dim>(i, j);
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        local_mat[loc_matrix_index( dim*n_shape_fct_+i , dim*n_shape_fct_+j )] += integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      { 
        // Integral_element - nabla phi_i \vec phi_j dx 
        // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
        integral = integrator.template integrate_vol_Dphiphi<hyEdge_dim>(i, j, dim);
        local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
        local_mat[loc_matrix_index(dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] -= integral;
    
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dim>(i, j, 2 * dim + 1);
        local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] += integral;
        // and from the penalty in the lower right diagonal block
        local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] 
          += tau * integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = integrator.template integrate_bdr_phiphi<hyEdge_dim>(i, j, 2 * dim + 0);
        local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
        // and from the penalty in the lower right diagonal block
        local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] 
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
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline std::array
< lSol_float_t, Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_ >
Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::assemble_rhs
(const std::array< std::array<lSol_float_t, n_shape_bdr_>, 2*hyEdge_dim >& lambda_values) const
{
  lSol_float_t integral;

  std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> right_hand_side;
  right_hand_side.fill(0.);

  hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
    hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
               "The size of lambda should be the amount of ansatz functions at boundary." );
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 0);
        right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
        right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j] * integral;
    
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 1);
        right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
        right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j] * integral;
      }
    }
  }
  
  return right_hand_side;
} // end of Diffusion_TensorialUniform::assemble_rhs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dim
>
Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::primal_at_boundary
( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
{
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 0);
        bdr_values[2*dim+0][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of Diffusion_TensorialUniform::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
inline std::array
< 
  std::array
  <
    lSol_float_t,
    Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dim
>
Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::dual_at_boundary
( const std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_>& coeffs ) const
{
  std::array<unsigned int, hyEdge_dim> dec_i;
  std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
  std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  lSol_float_t integral, integral1D;

  for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 0);
        bdr_values[2*dim+0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
        
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dim>(i, j, 2 * dim + 1);
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
< unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t >
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t >
std::array
<
  std::array
  <
    lSol_float_t,
    Hypercube<hyEdge_dim>::pow(sizeT)
  > ,
  Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()
>
Diffusion_TensorialUniform<hyEdge_dim,poly_deg,quad_deg,lSol_float_t>::bulk_values
(const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values) const
{
  std::array< lSol_float_t, n_loc_dofs_ > coefficients = solve_local_problem(lambda_values);

  std::array<std::array<lSol_float_t,Hypercube<hyEdge_dim>::pow(sizeT)>, system_dimension()> values;
  std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
  lSol_float_t fct_value;
 
  std::array<unsigned int, poly_deg+1> poly_indices;
  for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
  std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
    values1D = shape_fct_eval<lSol_float_t,Legendre>(poly_indices, abscissas);
      
  for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  { 
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
    {
      index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < system_dimension(); ++dim)
        values[dim][q] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
} // end of Diffusion_TensorialUniform::bulk_values


#endif // end of ifndef LSOL_DIFFUSION_HXX
