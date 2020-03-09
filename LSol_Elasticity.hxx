#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/ShapeFun1D.hxx>
#include <QuadratureTensorial.hxx>
#include <Hypercube.hxx>
#include <HyperHDG/LapackWrapper.hxx>

/*!*************************************************************************************************
 * \brief   Local solver for Poisson's equation on uniform hypergraph.
 *
 * \todo    Update doxygen in whole file!!!
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
< unsigned int hyEdge_dim, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t=double >
class ElasticRods_TensorialUniform
{
  public:
    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dim.
     * 
     * \retval  hyEdge_dim    Dimension of hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    /*!*********************************************************************************************
     * \brief   Decide whether gemetrical information is needed for local solver.
     * 
     * \retval  use_geom      True if geometrical information is used by local solver.
     **********************************************************************************************/
    static constexpr bool use_geometry() { return true; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number should be equal to \c n_dofs_per_nodeT of HyperNodeFactory.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return space_dim * Hypercube<hyEdge_dim-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return hyEdge_dim; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of quadrature points per spatial dimension.
     **********************************************************************************************/
    static constexpr unsigned int n_quads_1D_  = compute_n_quad_points<Gaussian>(quad_deg);
    /*!*********************************************************************************************
     * \brief   Number of local shape functions (with respect to all spatial dimensions).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_fct_ = Hypercube<hyEdge_dim>::pow(poly_deg + 1);
    /*!*********************************************************************************************
     * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     **********************************************************************************************/
    static constexpr unsigned int n_shape_bdr_ = Hypercube<hyEdge_dim-1>::pow(poly_deg + 1);
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
    assemble_loc_matrix ( const lSol_float_t tau )
    { 
      const std::array<lSol_float_t, n_quads_1D_> q_weights = quad_weights<quad_deg,Gaussian>();
      const std::array< std::array<lSol_float_t, n_quads_1D_ > , poly_deg + 1 > 
        trial(shape_fcts_at_quad_points<poly_deg, quad_deg, Gaussian, Legendre>()),
        deriv(shape_ders_at_quad_points<poly_deg, quad_deg, Gaussian, Legendre>());
      const std::array< std::array<lSol_float_t, 2> , poly_deg + 1 >
        trial_bdr(shape_fcts_at_bdrs<poly_deg, Legendre>());
  
      std::array<unsigned int, hyEdge_dim> dec_i, dec_j;
      lSol_float_t integral, integral1D;
      
      std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_> local_mat;
      local_mat.fill(0.);
  
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          index_decompose<hyEdge_dim>(j, poly_deg+1, dec_j);
          
          // Integral_element phi_i phi_j dx in diagonal blocks
          integral = 1.;
          for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_1D_; ++q)
              integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * trial[dec_i[dim_fct]][q];
            integral *= integral1D;
          }
          for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
            local_mat[loc_matrix_index( dim*n_shape_fct_+i , dim*n_shape_fct_+j )] += integral;
          
          for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          { 
            // Integral_element - nabla phi_i \vec phi_j dx 
            // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              integral1D = 0.;
              for (unsigned int q = 0; q < n_quads_1D_; ++q)
                integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * 
                  ( ( dim == dim_fct ) ? deriv[dec_i[dim_fct]][q] : trial[dec_i[dim_fct]][q] );
              integral *= integral1D;
            }
            local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
            local_mat[loc_matrix_index(dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] -= integral;
        
            // Corresponding boundary integrals from integration by parts in left lower blocks
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)
                integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] += integral;
            // Corresponding boundary integrals from integration by parts in left lower blocks
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)
                integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
        
            // Penalty in lower right diagonal block
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)
                integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] 
              += tau * integral;
            // Penalty in lower right diagonal block
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)
                integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+ j)]
              += tau * integral;
          }
        }
      }
  
      return local_mat;
    }
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
    
    /*!*********************************************************************************************
     * \brief  Do the pretprocessing to transfer global to local dofs.
     **********************************************************************************************/
    template<class GeomT> 
    inline std::array< std::array<double, n_shape_bdr_>, 2 * hyEdge_dim > node_dof_to_edge_dof
    ( const std::array< std::array<double, n_glob_dofs_per_node() >, 2 * hyEdge_dim > lambda,
      const GeomT& geom ) const
    {
      std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim > result;
      hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
      for (unsigned int i = 0; i < result.size(); ++i)
      {
        hy_assert( result[i].size() == 1 , "Only implemented in one dimension!" );
        result[i].fill(0.);
      }
  
      Point<space_dim> normal_vector = geom.normal(1);
  
      for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          result[i][0] += normal_vector[dim] * lambda[i][dim];
  
      return result;
    }
    /*!*********************************************************************************************
     * \brief  Do the postprocessing to transfer local to global dofs.
     **********************************************************************************************/
    template <class GeomT>
    inline std::array< std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dim >
    edge_dof_to_node_dof
    ( const std::array< std::array<double, n_shape_bdr_>, 2 * hyEdge_dim > lambda,
      const GeomT& geom ) const
    {
      hy_assert( n_shape_bdr_ == 1 , "This should be 1!")
      std::array< std::array<double, n_glob_dofs_per_node() > , 2*hyEdge_dim > result;
      Point<space_dim> normal_vector = geom.normal(1);
  
      for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          result[i][dim] = normal_vector[dim] * lambda[i][0];
  
      return result;
    }
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
    (const std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
    {
      std::array<unsigned int, hyEdge_dim> dec_i;
      std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      lSol_float_t integral, integral1D;
  
      std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> right_hand_side;
      right_hand_side.fill(0.);
  
      hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
                  "The size of the lambda values should be twice the dimension of a hyperedge." );
      for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
                   "The size of lambda should be the amount of ansatz functions at boundary." );
  
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
          for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          {
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
            right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j]*integral;
        
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
            right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j]*integral;
          }
        }
      }
  
      return right_hand_side;
    }
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    inline std::array<lSol_float_t, n_loc_dofs_ > solve_local_problem
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
    ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
    {
      std::array<unsigned int, hyEdge_dim> dec_i;
      std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
      lSol_float_t integral, integral1D;
    
      for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);
  
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
          for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          {
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+0][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
        
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+1][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
          }
        }
      }
      
      return bdr_values;
    }
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
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
            
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
          }
        }
      }
  
      return bdr_values;
    }
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
    ElasticRods_TensorialUniform(const constructor_value_type& tau)
    : tau_(tau), q_weights_(quad_weights<quad_deg, Gaussian>()),
      trial_(shape_fcts_at_quad_points<poly_deg, quad_deg, Gaussian, Legendre>()),
      trial_bdr_(shape_fcts_at_bdrs<poly_deg, Legendre>()),
      loc_mat_(assemble_loc_matrix(tau))
    { } 
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to matrix--vector multiplication.
     *
     * \todo    Include static asserts to check for correct dimensions.
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
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dim >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dim >&
      lambda_values, const GeomT& geom ) const
    {
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        lambda = node_dof_to_edge_dof(lambda_values, geom);
      std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda);
      
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > 
        bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      for (unsigned int i = 0; i < lambda.size(); ++i)
        for (unsigned int j = 0; j < lambda[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda[i][j];

      return edge_dof_to_node_dof(bdr_values, geom);
    }
    /*!*********************************************************************************************
     * \brief   Evaluate discrete function at given points.
     *
     * \todo    Include static asserts to check for correct dimensions.
     *
     * Function to evaluate primal variable of the solution at dyadic product of abscissas. The main
     * purpose of the function in closely related to plotting.
     *
     * \tparam  sizeT         Size of the passed \c std::array containing the abscissas.
     * \param   abscissas     Coordinates at whose tensor products the function is evaluated.
     * \param   lambda_values Coefficients of the associated skeletal function.
     * \retval  fct_val       Evaluation of dual variable at prescribed points.
     **********************************************************************************************/
    template<std::size_t sizeT, class GeomT>
    std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(sizeT)>
    primal_at_dyadic
    (const std::array<lSol_float_t, sizeT>& abscissas, const GeomT& geom,
     const std::array< std::array<lSol_float_t, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >&
      lambda_values) const
    {
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        lambda = node_dof_to_edge_dof(lambda_values, geom);
      std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda);
      
      std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(sizeT)> values;
      std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      lSol_float_t fct_value;

      std::array<unsigned int, poly_deg+1> poly_indices;
      for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
      std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
        values1D = shape_fct_eval<Legendre>(poly_indices, abscissas);
      
      values.fill(0.);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
        {
          index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
          fct_value = 1.;
          for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
          values[q] += coeffs[hyEdge_dim * n_shape_fct_ + i] * fct_value;
        }
      }

      return values;
    }
    /*!*********************************************************************************************
     * \brief   Evaluate discrete function at given points.
     *
     * \todo    Include static asserts to check for correct dimensions.
     *
     * Function to evaluate dual variable / part of the solution at dyadic product of abscissas. The
     * main purpose of the function in closely related to plotting.
     *
     * \tparam  sizeT         Size of the passed \c std::array containing the abscissas.
     * \param   abscissas     Coordinates at whose tensor products the function is evaluated.
     * \param   lambda_values Coefficients of the associated skeletal function.
     * \retval  fct_val       Evaluation of dual variable at prescribed points.
     **********************************************************************************************/
    template<std::size_t sizeT, class GeomT>
    std::array< std::array<lSol_float_t,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > 
    dual_at_dyadic
    (const std::array<lSol_float_t, sizeT>& abscissas, const GeomT& geom,
     const std::array< std::array<lSol_float_t, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >&
      lambda_values) const
    {
      std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        lambda = node_dof_to_edge_dof(lambda_values, geom);
      std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda);
      
      std::array< std::array<lSol_float_t, hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > values;
      std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      lSol_float_t fct_value;
      
      std::array<unsigned int, poly_deg+1> poly_indices;
      for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
      std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
        values1D = shape_fct_eval<Legendre>(poly_indices, abscissas);

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
          for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          values[q][dim] += coeffs[dim * n_shape_fct_ + i] * fct_value;
        }
      }
  
      return values;
    }
    
    
    template<typename AbscissaType, std::size_t AbscissaSize, class InputArrayType>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(AbscissaSize)>,system_dimension()>
    bulk_values (const std::array<AbscissaType,AbscissaSize>& abscissas,
	       const InputArrayType& lambda_values) const
    {
      std::array<std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(AbscissaSize)>,system_dimension()> result;
//      result[0] = primal_at_dyadic(abscissas, lambda_values, geom); ??
      return result;
    }
    
}; // end of class ElasticRods_TensorialUniform























//~ /*!*************************************************************************************************
 //~ * \brief   Local solver for Poisson's equation on uniform hypergraph.
 //~ *
 //~ * \todo    Update doxygen in whole file!!!
 //~ *
 //~ * This class contains the local solver for Poisson's equation, i.e.,
 //~ * \f[
 //~ *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 //~ * \f]
 //~ * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 //~ * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dim = 1) or hypergraph whose
 //~ * hyperedges are surfaces (\c hyEdge_dim = 2) or volumes (\c hyEdge_dim = 3). For this class, all
 //~ * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 //~ * information is needed by this class.
 //~ *
 //~ * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 //~ *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 //~ * \tparam  poly_deg      The polynomial degree of test and trial functions.
 //~ * \tparam  quad_deg      The order of the quadrature rule.
 //~ *
 //~ * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 //~ * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 //~ **************************************************************************************************/
//~ template
//~ < unsigned int hyEdge_dim, unsigned int space_dim, unsigned int poly_deg,
  //~ unsigned int quad_deg, typename lSol_float_t=double >
//~ class EulerBernoulliBeam
//~ {
  //~ public:
    //~ /*!*********************************************************************************************
     //~ * \brief   Return template parameter \c hyEdge_dim.
     //~ * 
     //~ * \retval  hyEdge_dim    Dimension of hypergraph's hyperedges.
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    //~ /*!*********************************************************************************************
     //~ * \brief   Decide whether gemetrical information is needed for local solver.
     //~ * 
     //~ * \retval  use_geom      True if geometrical information is used by local solver.
     //~ **********************************************************************************************/
    //~ static constexpr bool use_geometry() { return true; }
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate amount of global degrees of freedom per hypernode.
     //~ * 
     //~ * This number should be equal to \c n_dofs_per_nodeT of HyperNodeFactory.
     //~ *
     //~ * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int n_glob_dofs_per_node()
    //~ { return 2 * space_dim * Hypercube<hyEdge_dim-1>::pow(poly_deg + 1); }
    
    
    //~ static constexpr unsigned int node_value_dimension() { return 2*space_dim; }
    
    //~ static constexpr unsigned int system_dimension() { return 2*hyEdge_dim; }
    
    
  //~ private:
    //~ /*!*********************************************************************************************
     //~ * \brief   Number of quadrature points per spatial dimension.
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int n_quads_1D_  = compute_n_quad_points<Gaussian>(quad_deg);
    //~ /*!*********************************************************************************************
     //~ * \brief   Number of local shape functions (with respect to all spatial dimensions).
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int n_shape_fct_ = Hypercube<hyEdge_dim>::pow(poly_deg + 1);
    //~ /*!*********************************************************************************************
     //~ * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int n_shape_bdr_ = Hypercube<hyEdge_dim-1>::pow(poly_deg + 1);
    //~ /*!*********************************************************************************************
     //~ * \brief   Number of (local) degrees of freedom per hyperedge.
     //~ **********************************************************************************************/
    //~ static constexpr unsigned int n_loc_dofs_  = 2 * (hyEdge_dim+1) * n_shape_fct_;
    //~ /*!*********************************************************************************************
     //~ * \brief   Translate row and column indices to local index of entry in matrix.
     //~ * 
     //~ * Local \f$ n \times n \f$ matrices are encoded as arrays of size \f$n^2\f$. This function
     //~ * translates a row and a column index into the index of the long array, where the corresponding
     //~ * entry is located. Note that this is done column-wise (not row-wise as usually), to have the
     //~ * correct format for LAPACK.
     //~ *
     //~ * The function is static inline, since it is used in the constructor's initializer list.
     //~ *
     //~ * \param   row           Row index of local mtatrix entry.
     //~ * \param   column        Column index of local matrix entry.
     //~ * \retval  index         Overall index of local matrix entry.
     //~ **********************************************************************************************/
    //~ static inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column)
    //~ {
      //~ hy_assert( 0 <= row && row < n_loc_dofs_ ,
                 //~ "Row index must be >= 0 and smaller than total amount of local dofs." );
      //~ hy_assert( 0 <= column && column < n_loc_dofs_ ,
                 //~ "Column index must be >= 0 and smaller than total amount of local dofs." );

      //~ return column * n_loc_dofs_ + row;  // Transposed for LAPACK
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   Decompose index of local tensorial shape function into its (spatial) components.
     //~ *
     //~ * The function is static inline, since it is used in the constructor's initializer list.
     //~ *
     //~ * \param   index         Index of tensorial shape function.
     //~ * \param   range         Amount of 1D shape functions (the tensorial is made up from).
     //~ * \param   decomposition Array to be filled with the decomposition.
     //~ * \retval  decomposition Array filled with the decomposition.
     //~ **********************************************************************************************/
    //~ template<unsigned int dimT> static inline void index_decompose ( unsigned int index, 
      //~ unsigned int range, std::array<unsigned int, std::max(dimT,1U)>& decomposition )
    //~ {
      //~ if ( decomposition.size() == 1 )  decomposition[0] = index;
      //~ else
      //~ {
        //~ for (unsigned int dim = 0; dim < dimT; ++dim)
        //~ {
          //~ decomposition[dim] = index % range;
          //~ index /= range;
        //~ }
      //~ }
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief  Assemble local matrix for the local solver.
     //~ *
     //~ * The local solver neither depends on the geometry, nor on global functions. Thus, its local
     //~ * matrix is the same for all hyperedges and can be assembled once in the constructor. This is
     //~ * done in this function.
     //~ *
     //~ * The function is static inline, since it is used in the constructor's initializer list.
     //~ *
     //~ * \param   tau           Penalty parameter for HDG.
     //~ * \retval  loc_mat       Matrix of the local solver.
     //~ **********************************************************************************************/
    //~ static std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_ >
    //~ assemble_loc_matrix ( const lSol_float_t tau )
    //~ { 
      //~ const std::array<lSol_float_t, n_quads_1D_> q_weights = FuncQuad::quad_weights<quad_deg>();
      //~ const std::array< std::array<lSol_float_t, n_quads_1D_ > , poly_deg + 1 > 
        //~ trial(FuncQuad::shape_fcts_at_quad_points<poly_deg, quad_deg>()),
        //~ deriv(FuncQuad::shape_ders_at_quad_points<poly_deg, quad_deg>());
      //~ const std::array< std::array<lSol_float_t, 2> , poly_deg + 1 >
        //~ trial_bdr(FuncQuad::shape_fcts_at_bdrs<poly_deg>());
  
      //~ std::array<unsigned int, hyEdge_dim> dec_i, dec_j;
      //~ lSol_float_t integral, integral1D;
      
      //~ std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_> local_mat;
      //~ local_mat.fill(0.);
  
      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ {
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int j = 0; j < n_shape_fct_; ++j)
        //~ {
          //~ index_decompose<hyEdge_dim>(j, poly_deg+1, dec_j);
          
          //~ // Integral_element phi_i phi_j dx in diagonal blocks
          //~ integral = 1.;
          //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
          //~ {
            //~ integral1D = 0.;
            //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
              //~ integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * trial[dec_i[dim_fct]][q];
            //~ integral *= integral1D;
          //~ }
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ {
            //~ local_mat[loc_matrix_index( dim*n_shape_fct_+i , dim*n_shape_fct_+j )] += integral;
            //~ local_mat[loc_matrix_index( n_loc_dofs_/2 + dim*n_shape_fct_+i ,
                                        //~ n_loc_dofs_/2 + dim*n_shape_fct_+j )] += integral;
          //~ }
          //~ local_mat[loc_matrix_index( hyEdge_dim * n_shape_fct_ + i ,
                                      //~ n_loc_dofs_/2 + hyEdge_dim * n_shape_fct_ + j )] -= integral;
          
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ { 
            //~ // Integral_element - nabla phi_i \vec phi_j dx 
            //~ // = Integral_element - div \vec phi_i phi_j dx in right upper and left lower blocks
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ integral1D = 0.;
              //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                //~ integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * 
                  //~ ( ( dim == dim_fct ) ? deriv[dec_i[dim_fct]][q] : trial[dec_i[dim_fct]][q] );
              //~ integral *= integral1D;
            //~ }
            //~ local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
            //~ local_mat[loc_matrix_index(dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] -= integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + dim*n_shape_fct_+j)] -= integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+j)] -= integral;
        
            //~ // Corresponding boundary integrals from integration by parts in left lower blocks
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)
                //~ integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] += integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + dim*n_shape_fct_+j)] += integral;
            
            //~ // Corresponding boundary integrals from integration by parts in left lower blocks
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)
                //~ integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , dim*n_shape_fct_+j)] -= integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + dim*n_shape_fct_+j)] -= integral;
            
            //~ // Penalty in lower right diagonal block
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)
                //~ integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)] 
              //~ += tau * integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+j)] += tau * integral;
            
            //~ // Penalty in lower right diagonal block
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)
                //~ integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ local_mat[loc_matrix_index(hyEdge_dim*n_shape_fct_+i , hyEdge_dim*n_shape_fct_+j)]
              //~ += tau * integral;
            //~ local_mat[loc_matrix_index(n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+i ,
                                       //~ n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_+j)] += tau * integral;
          //~ }
        //~ }
      //~ }
  
      //~ return local_mat;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   (Globally constant) penalty parameter for HDG scheme.
     //~ **********************************************************************************************/
    //~ const lSol_float_t tau_;
    //~ /*!*********************************************************************************************
     //~ * \brief   Quadrature weights per spatial dimension.
     //~ **********************************************************************************************/
    //~ const std::array<lSol_float_t, n_quads_1D_> q_weights_;
    //~ /*!*********************************************************************************************
     //~ * \brief   Trial functions evaluated at quadrature points (per spatial dimensions).
     //~ **********************************************************************************************/
    //~ const std::array< std::array<lSol_float_t, n_quads_1D_ > , poly_deg + 1 > trial_;
    //~ /*!*********************************************************************************************
     //~ * \brief   Trial functions evaluated at boundaries {0,1} (per spatial dimension).
     //~ **********************************************************************************************/
    //~ const std::array< std::array<lSol_float_t, 2 > , poly_deg + 1 > trial_bdr_;
    //~ /*!*********************************************************************************************
     //~ * \brief   Local matrix for the local solver.
     //~ **********************************************************************************************/
    //~ const std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_ > loc_mat_;
    
    //~ /*!*********************************************************************************************
     //~ * \brief  Do the pretprocessing to transfer global to local dofs.
     //~ **********************************************************************************************/
    //~ template <class GeomT>
    //~ inline std::array< std::array<double, 2* n_shape_bdr_>, 2 * hyEdge_dim > node_dof_to_edge_dof
    //~ ( const std::array< std::array<double, n_glob_dofs_per_node() >, 2 * hyEdge_dim > lambda,
      //~ const GeomT& geom ) const
    //~ {
      //~ static_assert(space_dim == 2);
      //~ static_assert(hyEdge_dim == 1);
      //~ std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim > result;
      //~ hy_assert( result.size() == 4 , "Only implemented in one dimension!" );
      //~ for (unsigned int i = 0; i < result.size(); ++i)
      //~ {
        //~ hy_assert( result[i].size() == 1 , "Only implemented in one dimension!" );
        //~ result[i].fill(0.);
      //~ }
  
      //~ Point<space_dim> normal, normal_vector = geom.normal(1);
      //~ normal[0] = normal_vector[1]; normal[1] = -normal_vector[0];
  
      //~ for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        //~ for (unsigned int dim = 0; dim < space_dim; ++dim)
        //~ {
          //~ result[i][0] += normal[dim] * lambda[i][dim];
          //~ result[i+n_shape_bdr_][0] += normal[dim] * lambda[i+n_shape_bdr_][dim];
        //~ }
  
      //~ return result;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief  Do the postprocessing to transfer local to global dofs.
     //~ **********************************************************************************************/
    //~ template <class GeomT>
    //~ inline std::array< std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dim >
    //~ edge_dof_to_node_dof
    //~ ( const std::array< std::array<double, 2 * n_shape_bdr_>, 2 * hyEdge_dim > lambda,
      //~ const GeomT& geom ) const
    //~ {
      //~ hy_assert( n_shape_bdr_ == 1 , "This should be 1!")
      //~ std::array< std::array<double, n_glob_dofs_per_node() > , 2*hyEdge_dim > result;
      
      //~ Point<space_dim> normal, normal_vector = geom.normal(1);
      //~ normal[0] = normal_vector[1]; normal[1] = -normal_vector[0];
  
      //~ for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        //~ for (unsigned int dim = 0; dim < space_dim; ++dim)
        //~ {
          //~ result[i][dim] = normal[dim] * lambda[i][0];
          //~ result[i+n_shape_bdr_][dim] = normal[dim] * lambda[i+n_shape_bdr_][0];
        //~ }
  
      //~ return result;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief  Assemble local right hand for the local solver.
     //~ *
     //~ * The right hand side needs the values of the global degrees of freedom. Thus, it needs to be
     //~ * constructed individually for every hyperedge.
     //~ *
     //~ * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     //~ * \retval  loc_rhs       Local right hand side of the locasl solver.
     //~ **********************************************************************************************/
    //~ inline std::array<lSol_float_t, n_loc_dofs_ > assemble_rhs
    //~ (const std::array< std::array<lSol_float_t,2*n_shape_bdr_>, 2*hyEdge_dim >& lambda_values) const
    //~ {
      //~ std::array<unsigned int, hyEdge_dim> dec_i;
      //~ std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      //~ lSol_float_t integral, integral1D;
  
      //~ std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> right_hand_side;
      //~ right_hand_side.fill(0.);
  
      //~ hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
                  //~ "The size of the lambda values should be twice the dimension of a hyperedge." );
      //~ for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        //~ hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
                   //~ "The size of lambda should be the amount of ansatz functions at boundary." );
  
      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ {
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        //~ {
          //~ index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ {
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
            //~ right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j]*integral;
            //~ right_hand_side[n_loc_dofs_/2 + dim * n_shape_fct_ + i]
              //~ += lambda_values[2*dim+0][j+n_shape_bdr_] * integral;
            //~ right_hand_side[n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_ + i]
              //~ += tau_*lambda_values[2*dim+0][j+n_shape_bdr_]*integral;
        
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
            //~ right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j]*integral;
            //~ right_hand_side[n_loc_dofs_/2 + dim*n_shape_fct_ + i]
              //~ -= lambda_values[2*dim+1][j+n_shape_bdr_] * integral;
            //~ right_hand_side[n_loc_dofs_/2 + hyEdge_dim*n_shape_fct_ + i]
              //~ += tau_*lambda_values[2*dim+1][j+n_shape_bdr_]*integral;
          //~ }
        //~ }
      //~ }
  
      //~ return right_hand_side;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief  Solve local problem.
     //~ *
     //~ * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     //~ * \retval  loc_sol       Solution of the local problem.
     //~ **********************************************************************************************/
    //~ inline std::array<lSol_float_t, n_loc_dofs_ > solve_local_problem
    //~ (const std::array< std::array<lSol_float_t,2*n_shape_bdr_>, 2*hyEdge_dim >& lambda_values) const
    //~ {
      //~ // A copy of loc_mat_ is created, since LAPACK will destroy the matrix values.
      //~ std::array<lSol_float_t, n_loc_dofs_ * n_loc_dofs_> local_matrix = loc_mat_;
      //~ // The local right hand side is assembled (and will also be destroyed by LAPACK).
      //~ std::array<lSol_float_t, n_loc_dofs_> right_hand_side = assemble_rhs(lambda_values);
      //~ // LAPACK solves local_matix * return_value = right_hand_side.
      //~ try { return lapack_solve<n_loc_dofs_>(local_matrix, right_hand_side); }
      //~ catch (LASolveException& exc)
      //~ {
        //~ hy_assert( 0 == 1 ,
                   //~ exc.what() << std::endl << "This can happen if quadrature is too inaccurate!" );
        //~ throw exc;
      //~ }
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate primal variable at boundary.
     //~ *
     //~ * Function to evaluate primal variable of the solution. This function is needed to calculate
     //~ * the local numerical fluxes.
     //~ *
     //~ * \param   coeffs        Coefficients of the local solution.
     //~ * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     //~ **********************************************************************************************/
    //~ inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > primal_at_boundary
    //~ ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
    //~ {
      //~ std::array<unsigned int, hyEdge_dim> dec_i;
      //~ std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
      //~ lSol_float_t integral, integral1D;
    
      //~ for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);
  
      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ {
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        //~ {
          //~ index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ {
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ bdr_values[2*dim+0][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
            //~ bdr_values[2*dim+0][j+n_shape_bdr_]
              //~ += coeffs[n_loc_dofs_/2 + hyEdge_dim * n_shape_fct_ + i] * integral;
        
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ bdr_values[2*dim+1][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
            //~ bdr_values[2*dim+1][j+n_shape_bdr_]
              //~ += coeffs[n_loc_dofs_/2 + hyEdge_dim * n_shape_fct_ + i] * integral;
          //~ }
        //~ }
      //~ }
      
      //~ return bdr_values;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate dual variable at boundary.
     //~ *
     //~ * Function to evaluate dual variable of the solution. This function is needed to calculate the
     //~ * local numerical fluxes.
     //~ *
     //~ * \param   coeffs        Coefficients of the local solution.
     //~ * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     //~ **********************************************************************************************/
    //~ inline std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > dual_at_boundary
    //~ ( const std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_>& coeffs ) const
    //~ {
      //~ std::array<unsigned int, hyEdge_dim> dec_i;
      //~ std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
      //~ lSol_float_t integral, integral1D;
    
      //~ for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);

      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ {
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        //~ {
          //~ index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ {
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ bdr_values[2*dim+0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
            //~ bdr_values[2*dim+0][j+n_shape_bdr_]
              //~ -= coeffs[n_loc_dofs_/2 + dim * n_shape_fct_ + i] * integral;
            
            //~ integral = 1.;
            //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ {
              //~ if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              //~ else
              //~ {
                //~ integral1D = 0.;
                //~ for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  //~ integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              //~ }
              //~ integral *= integral1D;
            //~ }
            //~ bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
            //~ bdr_values[2*dim+1][j+n_shape_bdr_]
              //~ += coeffs[n_loc_dofs_/2 + dim * n_shape_fct_ + i] * integral;
          //~ }
        //~ }
      //~ }
  
      //~ return bdr_values;
    //~ }
  //~ public:
    //~ /*!*********************************************************************************************
     //~ * \brief   Class is constructed using a single double indicating the penalty parameter.
     //~ **********************************************************************************************/
    //~ typedef lSol_float_t constructor_value_type;
    //~ /*!*********************************************************************************************
     //~ * \brief   Constructor for local solver.
     //~ *
     //~ * \param   tau           Penalty parameter of HDG scheme.
     //~ **********************************************************************************************/
    //~ EulerBernoulliBeam(const constructor_value_type& tau)
    //~ : tau_(tau), q_weights_(FuncQuad::quad_weights<quad_deg>()),
      //~ trial_(FuncQuad::shape_fcts_at_quad_points<poly_deg, quad_deg>()),
      //~ trial_bdr_(FuncQuad::shape_fcts_at_bdrs<poly_deg>()),
      //~ loc_mat_(assemble_loc_matrix(tau))
    //~ { } 
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate local contribution to matrix--vector multiplication.
     //~ *
     //~ * \todo    Include static asserts to check for correct dimensions.
     //~ *
     //~ * Execute matrix--vector multiplication y = A * x, where x represents the vector containing the
     //~ * skeletal variable (adjacent to one hyperedge), and A is the condensed matrix arising from the
     //~ * HDG discretization. This function does this multiplication (locally) for one hyperedge. The
     //~ * hyperedge is no parameter, since all hyperedges are assumed to have the same properties.
     //~ *
     //~ * \param   lambda_values Local part of vector x.
     //~ * \retval  vecAx         Local part of vector A * x.
     //~ **********************************************************************************************/
    //~ template <class GeomT>
    //~ std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dim >
    //~ numerical_flux_from_lambda
    //~ (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dim >&
      //~ lambda_values, const GeomT& geom ) const
    //~ {
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        //~ lambda = node_dof_to_edge_dof(lambda_values, geom);
      //~ std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda);
      
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2 * hyEdge_dim > 
        //~ bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      //~ for (unsigned int i = 0; i < lambda.size(); ++i)
        //~ for (unsigned int j = 0; j < lambda[i].size(); ++j)
          //~ bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda[i][j];

      //~ return edge_dof_to_node_dof(bdr_values, geom);
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate discrete function at given points.
     //~ *
     //~ * \todo    Include static asserts to check for correct dimensions.
     //~ *
     //~ * Function to evaluate primal variable of the solution at dyadic product of abscissas. The main
     //~ * purpose of the function in closely related to plotting.
     //~ *
     //~ * \tparam  sizeT         Size of the passed \c std::array containing the abscissas.
     //~ * \param   abscissas     Coordinates at whose tensor products the function is evaluated.
     //~ * \param   lambda_values Coefficients of the associated skeletal function.
     //~ * \retval  fct_val       Evaluation of dual variable at prescribed points.
     //~ **********************************************************************************************/
    //~ template<std::size_t sizeT, class GeomT>
    //~ std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(sizeT)>
    //~ primal_at_dyadic
    //~ (const std::array<lSol_float_t, sizeT>& abscissas, const GeomT& geom,
     //~ const std::array< std::array<lSol_float_t, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >&
      //~ lambda_values) const
    //~ {
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        //~ lambda = node_dof_to_edge_dof(lambda_values, geom);
      //~ std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda);
      
      //~ std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(sizeT)> values;
      //~ std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      //~ lSol_float_t fct_value;

      //~ std::array<unsigned int, poly_deg+1> poly_indices;
      //~ for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
      //~ std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
        //~ values1D = FuncQuad::shape_fct_eval(poly_indices, abscissas);
      
      //~ values.fill(0.);
      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ {
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
        //~ {
          //~ index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
          //~ fct_value = 1.;
          //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
          //~ values[q] += coeffs[hyEdge_dim * n_shape_fct_ + i] * fct_value;
        //~ }
      //~ }

      //~ return values;
    //~ }
    //~ /*!*********************************************************************************************
     //~ * \brief   Evaluate discrete function at given points.
     //~ *
     //~ * \todo    Include static asserts to check for correct dimensions.
     //~ *
     //~ * Function to evaluate dual variable / part of the solution at dyadic product of abscissas. The
     //~ * main purpose of the function in closely related to plotting.
     //~ *
     //~ * \tparam  sizeT         Size of the passed \c std::array containing the abscissas.
     //~ * \param   abscissas     Coordinates at whose tensor products the function is evaluated.
     //~ * \param   lambda_values Coefficients of the associated skeletal function.
     //~ * \retval  fct_val       Evaluation of dual variable at prescribed points.
     //~ **********************************************************************************************/
    //~ template<std::size_t sizeT, class GeomT>
    //~ std::array< std::array<lSol_float_t,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > 
    //~ dual_at_dyadic
    //~ (const std::array<lSol_float_t, sizeT>& abscissas, const GeomT& geom,
     //~ const std::array< std::array<lSol_float_t, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >&
      //~ lambda_values) const
    //~ {
      //~ std::array< std::array<lSol_float_t, n_shape_bdr_> , 2*hyEdge_dim >
        //~ lambda = node_dof_to_edge_dof(lambda_values, geom);
      //~ std::array<lSol_float_t, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda);
      
      //~ std::array< std::array<lSol_float_t, hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > values;
      //~ std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      //~ lSol_float_t fct_value;
      
      //~ std::array<unsigned int, poly_deg+1> poly_indices;
      //~ for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
      //~ std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
        //~ values1D = FuncQuad::shape_fct_eval(poly_indices, abscissas);

      //~ for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
      //~ for (unsigned int i = 0; i < n_shape_fct_; ++i)
      //~ { 
        //~ index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
        //~ for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
        //~ {
          //~ index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
          //~ fct_value = 1.;
          //~ for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            //~ fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
          //~ for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
          //~ values[q][dim] += coeffs[dim * n_shape_fct_ + i] * fct_value;
        //~ }
      //~ }
  
      //~ return values;
    //~ }
    
    
    //~ template<typename AbscissaType, std::size_t AbscissaSize, class InputArrayType>
    //~ std::array<std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(AbscissaSize)>,system_dimension()>
    //~ bulk_values (const std::array<AbscissaType,AbscissaSize>& abscissas,
	       //~ const InputArrayType& lambda_values) const
    //~ {
      //~ std::array<std::array<lSol_float_t, Hypercube<hyEdge_dim>::pow(AbscissaSize)>,system_dimension()> result;
//~ //      result[0] = primal_at_dyadic(abscissas, lambda_values, geom); ??
      //~ return result;
    //~ }
    
//~ }; // end of class ElasticRods_TensorialUniform
