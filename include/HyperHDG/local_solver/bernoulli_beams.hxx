#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/shape_fun_1d.hxx>
#include <HyperHDG/quadrature_tensorial.hxx>

#include <HyperHDG/local_solver/diffusion.hxx>
#include <HyperHDG/local_solver/bilaplacian.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/tensorial_shape_fun.hxx>

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
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t = double,
  typename diffusion_sol_t = DiffusionUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>
>
class LengtheningBeam
{
  public:

    typedef struct empty_class {} data_type;

    typedef lSol_float_t solver_float_t;

    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dimT.
     * 
     * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
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
    { return 2 * space_dim * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return space_dim; }
    
  private:
    const diffusion_sol_t diffusion;
    /*!*********************************************************************************************
     * \brief  Do the pretprocessing to transfer global to local dofs.
     **********************************************************************************************/
    template<class hyEdgeT> 
    inline std::array< std::array<double,diffusion_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT >
    node_dof_to_edge_dof
    ( const std::array< std::array<double, n_glob_dofs_per_node() >, 2 * hyEdge_dimT >& lambda,
      hyEdgeT& hyper_edge ) const
    {
      std::array<std::array<double, diffusion_sol_t::n_glob_dofs_per_node()>,2*hyEdge_dimT> result;
      hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
      for (unsigned int i = 0; i < result.size(); ++i)
      {
        hy_assert( result[i].size() == 1 , "Only implemented in one dimension!" );
        result[i].fill(0.);
      }
  
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.inner_normal(1);
  
      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          result[i][0] += normal_vector[dim] * lambda[i][dim];
  
      return result;
    }
    /*!*********************************************************************************************
     * \brief  Do the postprocessing to transfer local to global dofs.
     **********************************************************************************************/
    template <class hyEdgeT>
    inline std::array< std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT >
    edge_dof_to_node_dof
    ( const std::array<std::array<double,diffusion_sol_t::n_glob_dofs_per_node()>,2*hyEdge_dimT>& 
      lambda, hyEdgeT& hyper_edge ) const
    {
      hy_assert( diffusion_sol_t::n_glob_dofs_per_node() == 1 , "This should be 1!" );
      std::array< std::array<double, n_glob_dofs_per_node() > , 2*hyEdge_dimT > result;
      for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.inner_normal(1);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          result[i][dim] = normal_vector[dim] * lambda[i][0];
  
      return result;
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
    LengtheningBeam(const constructor_value_type& tau = 1.) : diffusion(tau)  { }
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
    template <class hyEdgeT>
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dimT >&
      lambda_values, hyEdgeT& hyper_edge, const lSol_float_t time = 0.
    ) const
    {
      std::array<std::array<lSol_float_t,diffusion_sol_t::n_glob_dofs_per_node()>,2*hyEdge_dimT >
        lambda = node_dof_to_edge_dof(lambda_values, hyper_edge);

      return edge_dof_to_node_dof(diffusion.numerical_flux_from_lambda(lambda),hyper_edge);
    }

    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      LengtheningBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()
    >
    bulk_values
    (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
     hyEdgeT& hyper_edge, const lSol_float_t time = 0.) const
    {
      std::array
      <
        std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
        LengtheningBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()
      > result;

      auto bulk = diffusion.bulk_values(abscissas,node_dof_to_edge_dof(lambda_values, hyper_edge));
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.inner_normal(1);

      for(unsigned int dim = 0; dim < result.size(); ++dim)
        for (unsigned int q = 0; q < result[dim].size(); ++q)
          result[dim][q] = bulk[1][q] * normal_vector[dim];

      return result;
    }

  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1> result;
    return result;
  }
}; // end of class LengtheningBeam


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
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t = double,
  typename bilaplacian_sol_t = BilaplacianUniform<hyEdge_dimT,poly_deg,quad_deg,lSol_float_t>
>
class BernoulliBendingBeam
{
  public:
    typedef struct empty_class {} data_type;
    typedef lSol_float_t solver_float_t;

    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dimT.
     * 
     * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Specify whether advanced functuions are implemented for this class
     * 
     * Advanced functions are numerical_flux_from_rhs, dirichlet_coeffs, and neumann_coeffs.
     * 
     * \todo    Replace this by concept in C++20!
     **********************************************************************************************/    
    static constexpr bool advanced_functions() { return false; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * space_dim * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return space_dim; }
    
    
  private:
    const bilaplacian_sol_t bilaplacian_solver;   
    /*!*********************************************************************************************
     * \brief  Do the pretprocessing to transfer global to local dofs.
     **********************************************************************************************/
    template<class hyEdgeT> 
    inline std::array< std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT > node_dof_to_edge_dof
    ( const std::array< std::array<double, n_glob_dofs_per_node() >, 2 * hyEdge_dimT >& lambda,
      hyEdgeT& hyper_edge, const unsigned int outer_index ) const
    {
      std::array< std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()> , 2*hyEdge_dimT > result;
      hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
      for (unsigned int i = 0; i < result.size(); ++i)
      {
        hy_assert( result[i].size() == 2 , "Only implemented in one dimension!" );
        result[i].fill(0.);
      }
  
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(outer_index);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          result[i][0] += normal_vector[dim] * lambda[i][dim];
          result[i][1] += normal_vector[dim] * lambda[i][space_dim + dim];
        }
  
      return result;
    }
    /*!*********************************************************************************************
     * \brief  Do the postprocessing to transfer local to global dofs.
     **********************************************************************************************/
    template <class hyEdgeT>
    inline std::array< std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT >
    edge_dof_to_node_dof
    ( const std::array< std::array<double, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2 * hyEdge_dimT >& lambda,
      hyEdgeT& hyper_edge, const unsigned int outer_index ) const
    {
      hy_assert( bilaplacian_sol_t::n_glob_dofs_per_node() == 2 , "This should be 1*2!")
      std::array< std::array<double, n_glob_dofs_per_node() > , 2*hyEdge_dimT > result;
      for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(outer_index);
  
      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          result[i][dim] = normal_vector[dim] * lambda[i][0];
          result[i][space_dim + dim] = normal_vector[dim] * lambda[i][1];
        }
  
      return result;
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
    BernoulliBendingBeam(const constructor_value_type& tau = 1.) : bilaplacian_solver(tau)  { }
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
    template <class hyEdgeT>
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dimT >&
      lambda_values, hyEdgeT& hyper_edge, const lSol_float_t time = 0. ) const
    {
      std::array<std::array<lSol_float_t,bilaplacian_sol_t::n_glob_dofs_per_node()>,2*hyEdge_dimT>
        lambda;

      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > result, aux;
      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)  result[i].fill(0.);

      for (unsigned int dim = 0; dim < space_dim - hyEdge_dimT; ++dim)
      {
        lambda = node_dof_to_edge_dof(lambda_values, hyper_edge, dim);
        aux = edge_dof_to_node_dof(bilaplacian_solver.numerical_flux_from_lambda(lambda), hyper_edge, dim);

        for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
          for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
            result[i][j] += aux[i][j];
      }

      return result;
    }
    
    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to right-hand side vector by global right-hand side.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()>, 2*hyEdge_dimT > numerical_flux_from_rhs
    ( const std::array< unsigned int, 2*hyEdge_dimT > & bound_cond, hyEdgeT& hyper_edge )  const
    {
      std::array< std::array<lSol_float_t, bilaplacian_sol_t::n_glob_dofs_per_node()> , 2 * hyEdge_dimT > bdr_values;            
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      BernoulliBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()
    >
    bulk_values
    ( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
      hyEdgeT& hyper_edge, const lSol_float_t time = 0. ) const
    {
      std::array<std::array<lSol_float_t,Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()> values;
      for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
      for (unsigned int dim_on = 0; dim_on < space_dim - hyEdge_dimT; ++dim_on)
      {
        auto bulk = bilaplacian_solver.bulk_values(abscissas,node_dof_to_edge_dof(lambda_values, hyper_edge,dim_on));
        Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(dim_on);

        for(unsigned int dim = 0; dim < values.size(); ++dim)
          for (unsigned int q = 0; q < values[dim].size(); ++q)
            values[dim][q] += bulk[1][q] * normal_vector[dim];
      }
  
      return values;
    }

  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1> result;
    return result;
  }

}; // end of class BernoulliBendingBeam


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
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t = double
>
class LengtheningBernoulliBendingBeam
{
  public:
    typedef struct empty_class {} data_type;
    typedef lSol_float_t solver_float_t;
    
    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dimT.
     * 
     * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Specify whether advanced functuions are implemented for this class
     * 
     * Advanced functions are numerical_flux_from_rhs, dirichlet_coeffs, and neumann_coeffs.
     * 
     * \todo    Replace this by concept in C++20!
     **********************************************************************************************/    
    static constexpr bool advanced_functions() { return false; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * space_dim * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return space_dim; }
    
    
  private:
    
    const lSol_float_t tau_;

    const LengtheningBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t> len_beam;
    const BernoulliBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t> ben_beam;

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
    LengtheningBernoulliBendingBeam(const constructor_value_type& tau = 1.)
    : tau_(tau), len_beam(tau), ben_beam(tau)  { }
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
    template <class hyEdgeT>
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dimT >&
      lambda_values, hyEdgeT& hyper_edge, const lSol_float_t time = 0. ) const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > result, aux;
      
      result = len_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);
      aux = ben_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
          result[i][j] += aux[i][j];

      return result;
    }

    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to right-hand side vector by global right-hand side.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()>, 2*hyEdge_dimT > numerical_flux_from_rhs
    ( const std::array< unsigned int, 2*hyEdge_dimT > & bound_cond, hyEdgeT& hyper_edge, const lSol_float_t = 0. )  const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > bdr_values;            
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()>
    bulk_values
    (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
     hyEdgeT& hyper_edge, const lSol_float_t time = 0.) const
    {
      std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()>
        result, auxiliary;
      result = len_beam.bulk_values(abscissas,lambda_values,hyper_edge);
      auxiliary = ben_beam.bulk_values(abscissas,lambda_values,hyper_edge);

      for (unsigned int i = 0; i < system_dimension(); ++i)
        for (unsigned int j = 0; j < Hypercube<hyEdge_dimT>::pow(sizeT); ++j)
          result[i][j] += auxiliary[i][j];

      return result;
    }

  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1> result;
    for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
    return result;
  }

}; // end of class LengtheningBernoulliBendingBeam










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
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
<
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t = double
>
class TimoschenkoBendingBeam
{
  public:
    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dimT.
     * 
     * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Specify whether advanced functuions are implemented for this class
     * 
     * Advanced functions are numerical_flux_from_rhs, dirichlet_coeffs, and neumann_coeffs.
     * 
     * \todo    Replace this by concept in C++20!
     **********************************************************************************************/    
    static constexpr bool advanced_functions() { return false; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 3 * space_dim * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return space_dim; }
    
    
  private:
    /*!*********************************************************************************************
     * \brief   Number of quadrature points per spatial dimension.
     **********************************************************************************************/
    static constexpr unsigned int n_quads_1D_  = compute_n_quad_points<Gaussian>(quad_deg);
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
    static constexpr unsigned int n_loc_dofs_  = (2*hyEdge_dimT+2) * n_shape_fct_;
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
     * \brief  Do the pretprocessing to transfer global to local dofs.
     **********************************************************************************************/
    template<class hyEdgeT> 
    inline std::array< std::array<double, 3 * n_shape_bdr_>, 2 * hyEdge_dimT > node_dof_to_edge_dof
    ( const std::array< std::array<double, n_glob_dofs_per_node() >, 2 * hyEdge_dimT >& lambda,
      hyEdgeT& hyper_edge, const unsigned int outer_index ) const
    {
      std::array< std::array<double, 2 * n_shape_bdr_> , 2*hyEdge_dimT > result;
      hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
      for (unsigned int i = 0; i < result.size(); ++i)
      {
        hy_assert( result[i].size() == 2 , "Only implemented in one dimension!" );
        result[i].fill(0.);
      }
  
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(outer_index);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          result[i][0] += normal_vector[dim] * lambda[i][0 * space_dim + dim];  // Displacement
          result[i][1] += normal_vector[dim] * lambda[i][1 * space_dim + dim];  // Laplacian
          result[i][2] += normal_vector[dim] * lambda[i][2 * space_dim + dim];  // Torsion
        }
  
      return result;
    }
    /*!*********************************************************************************************
     * \brief  Do the postprocessing to transfer local to global dofs.
     **********************************************************************************************/
    template <class hyEdgeT>
    inline std::array< std::array<double, n_glob_dofs_per_node()>, 2 * hyEdge_dimT >
    edge_dof_to_node_dof
    ( const std::array< std::array<double, 2 * n_shape_bdr_>, 2 * hyEdge_dimT >& lambda,
      hyEdgeT& hyper_edge, const unsigned int outer_index ) const
    {
      hy_assert( n_shape_bdr_ == 1 , "This should be 1!")
      std::array< std::array<double, n_glob_dofs_per_node() > , 2*hyEdge_dimT > result;
      for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
      Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(outer_index);
  
      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          result[i][0 * space_dim + dim] = normal_vector[dim] * lambda[i][0];  // Displacement
          result[i][1 * space_dim + dim] = normal_vector[dim] * lambda[i][1];  // Laplacian
          result[i][2 * space_dim + dim] = normal_vector[dim] * lambda[i][2];  // Torsion
        }
  
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
    inline SmallVec< n_loc_dofs_, lSol_float_t > assemble_rhs
    ( const std::array<std::array<lSol_float_t, 3 * n_shape_bdr_>, 2*hyEdge_dimT>& lambda_values )
    const;
    
    /*!*********************************************************************************************
     * \brief  Solve local problem.
     *
     * \param   lambda_values Global degrees of freedom associated to the hyperedge.
     * \retval  loc_sol       Solution of the local problem.
     **********************************************************************************************/
    inline std::array< lSol_float_t, n_loc_dofs_ > solve_local_problem
    ( const std::array<std::array<lSol_float_t, 3 * n_shape_bdr_> , 2*hyEdge_dimT>& lambda_values )
    const
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
    inline std::array< std::array<lSol_float_t, 3 * n_shape_bdr_> , 2 * hyEdge_dimT >
    primal_at_boundary ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const;
    /*!*********************************************************************************************
     * \brief   Evaluate dual variable at boundary.
     *
     * Function to evaluate dual variable of the solution. This function is needed to calculate the
     * local numerical fluxes.
     *
     * \param   coeffs        Coefficients of the local solution.
     * \retval  bdr_coeffs    Coefficients of respective (dim-1) dimensional function at boundaries.
     **********************************************************************************************/
    inline std::array< std::array<lSol_float_t, 3 * n_shape_bdr_> , 2 * hyEdge_dimT >
    dual_at_boundary ( const std::array<lSol_float_t, n_loc_dofs_>& coeffs ) const;
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
    TimoschenkoBendingBeam(const constructor_value_type& tau = 1.)
    : tau_(tau), loc_mat_(assemble_loc_matrix(tau))  { }
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
    template <class hyEdgeT>
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dimT >&
      lambda_values, hyEdgeT& hyper_edge ) const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > result, aux;
      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)  result[i].fill(0.);

      for (unsigned int dim = 0; dim < space_dim - hyEdge_dimT; ++dim)
      {
        std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2*hyEdge_dimT >
          lambda = node_dof_to_edge_dof(lambda_values, hyper_edge, dim);
        std::array<lSol_float_t, n_loc_dofs_ > coeffs = solve_local_problem(lambda);
      
        std::array< std::array<lSol_float_t, 2*n_shape_bdr_> , 2 * hyEdge_dimT > 
          bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
        for (unsigned int i = 0; i < lambda.size(); ++i)
          for (unsigned int j = 0; j < lambda[i].size(); ++j)
            bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda[i][j];

        aux = edge_dof_to_node_dof(bdr_values, hyper_edge, dim);
        for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
          for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
            result[i][j] += aux[i][j];
      }

      return result;
    }

    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to right-hand side vector by global right-hand side.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()>, 2*hyEdge_dimT > numerical_flux_from_rhs
    ( const std::array< unsigned int, 2*hyEdge_dimT > & bound_cond, hyEdgeT& hyper_edge )  const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > bdr_values;            
      return bdr_values;
    }
    
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
    std::array
    <
      std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,
      TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::
        system_dimension()
    >
    bulk_values
    ( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
      hyEdgeT& hyper_edge, const lSol_float_t time = 0. ) const;

  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const;

}; // end of class TimoschenkoBendingBeam



// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF TimoschenkoBendingBeam
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix: With artificial $\Delta u = 0 \text{ on } \partial \elem$
// -------------------------------------------------------------------------------------------------

template
<
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t
>
SmallSquareMat
<
  TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_,
  lSol_float_t
>
TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::
assemble_loc_matrix ( const lSol_float_t tau )
{ 
  constexpr unsigned int n_dofs_lap = (hyEdge_dimT+3) * n_shape_fct_;
  const IntegratorTensorial<poly_deg,quad_deg,Gaussian,Legendre,lSol_float_t> integrator;
  lSol_float_t integral;
  
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = integrator.template integrate_vol_phiphi<hyEdge_dimT>(i, j);
      local_mat(              hyEdge_dimT * n_shape_fct_ + i ,
                 n_dofs_lap                              + j ) += integral;
      local_mat( n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i ,
                                                         + j ) -= integral;
      local_mat( n_dofs_lap + hyEdge_dimT * n_shape_fct_ + i ,
                 n_dofs_lap + hyEdge_dimT * n_shape_fct_ + j ) += integral;
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
} // end of TimoschenkoBendingBeam::assemble_loc_matrix


// -------------------------------------------------------------------------------------------------
// assemble_rhs
// -------------------------------------------------------------------------------------------------

template
<
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t
>
inline SmallVec
< 
  TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::n_loc_dofs_,
  lSol_float_t
>
TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::assemble_rhs
(const std::array< std::array<lSol_float_t, 3 * n_shape_bdr_>, 2*hyEdge_dimT >& lambda_values) const
{
  constexpr unsigned int n_dofs_lap = (hyEdge_dimT+1) * n_shape_fct_;
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
        
        right_hand_side[n_dofs_lap + i] -= lambda_values[2*dim+0][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + n_shape_fct_ + i]
          += lambda_values[2*dim+0][2 * n_shape_bdr_ + j] * integral;
    
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
        right_hand_side[hyEdge_dimT*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j] * integral;
        
        right_hand_side[n_dofs_lap + i] += lambda_values[2*dim+1][n_shape_bdr_ + j] * integral;
        right_hand_side[n_dofs_lap + n_shape_fct_ + i]
          -= lambda_values[2*dim+1][2 * n_shape_bdr_ + j] * integral;
      }
    }
  }
  
  return right_hand_side;
} // end of TimoschenkoBendingBeam::assemble_rhs


// -------------------------------------------------------------------------------------------------
// primal_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t
>
inline std::array
< 
  std::array
  <
    lSol_float_t,
    3 * TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::
primal_at_boundary ( const std::array<lSol_float_t, n_loc_dofs_ >& coeffs ) const
{
  constexpr unsigned int n_dofs_lap = (hyEdge_dimT+1) * n_shape_fct_;
  std::array< std::array<lSol_float_t, 3 * n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
  lSol_float_t integral;

  for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dimT; ++dim_n)  bdr_values[dim_n].fill(0.);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 0);
        bdr_values[2*dim+0][0*n_shape_bdr_ + j] += coeffs[hyEdge_dimT*n_shape_fct_ + i] * integral;
        bdr_values[2*dim+0][1*n_shape_bdr_ + j] += coeffs[n_dofs_lap + i] * integral;
        bdr_values[2*dim+0][2*n_shape_bdr_ + j] += coeffs[n_dofs_lap + n_shape_bdr_ + i] * integral;

        integral = integrator.template integrate_bdr_phipsi<hyEdge_dimT>(i, j, 2 * dim + 1);
        bdr_values[2*dim+1][0*n_shape_bdr_ + j] += coeffs[hyEdge_dimT*n_shape_fct_ + i] * integral;
        bdr_values[2*dim+1][1*n_shape_bdr_ + j] += coeffs[n_dofs_lap + i] * integral;
        bdr_values[2*dim+1][2*n_shape_bdr_ + j] += coeffs[n_dofs_lap + n_shape_bdr_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
} // end of TimoschenkoBendingBeam::primal_at_boundary


// -------------------------------------------------------------------------------------------------
// dual_at_boundary
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t
>
inline std::array
< 
  std::array
  <
    lSol_float_t,
    3 * TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::n_shape_bdr_
  > ,
  2 * hyEdge_dimT
>
TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::
dual_at_boundary ( const std::array<lSol_float_t, n_loc_dofs_>& coeffs ) const
{
  constexpr unsigned int n_dofs_lap = (hyEdge_dimT+1) * n_shape_fct_;
  std::array< std::array<lSol_float_t, 3 * n_shape_bdr_> , 2 * hyEdge_dimT > bdr_values;
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
} // end of TimoschenkoBendingBeam::dual_at_boundary


// -------------------------------------------------------------------------------------------------
// bulk_values
// -------------------------------------------------------------------------------------------------

template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg,
  typename lSol_float_t
>
template < typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT >
std::array
<
  std::array
  <
    lSol_float_t,
    Hypercube<hyEdge_dimT>::pow(sizeT)
  > ,
  TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::system_dimension()
>
TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t>::bulk_values
( const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
  hyEdgeT& hyper_edge, const lSol_float_t time ) const
{
  std::array<std::array<lSol_float_t,Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()> values;
  for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
  for (unsigned int dim_on = 0; dim_on < space_dim - hyEdge_dimT; ++dim_on)
  {
    std::array< lSol_float_t, n_loc_dofs_ > coefficients
      = solve_local_problem(node_dof_to_edge_dof(lambda_values, hyper_edge, dim_on));

    std::array<std::array<lSol_float_t,Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()> values
        = sum_all_in_tensorial_points_with_coefficients<lSol_float_t,abscissa_float_t,Legendre,
                                                        hyEdge_dimT,abscissas.size(), poly_deg, system_dimension()>(abscissas, coefficients);
    std::array<unsigned int, hyEdge_dimT> dec_i, dec_q;
    lSol_float_t fct_value;
 
    std::array<unsigned int, poly_deg+1> poly_indices;
    for (unsigned int i = 0; i < poly_deg+1; ++i) poly_indices[i] = i;
    Point<space_dim,lSol_float_t> normal_vector = hyper_edge.geometry.outer_normal(dim_on);
    std::array< std::array<lSol_float_t, abscissas.size()>, poly_deg+1 > 
    values1D = shape_fct_eval<lSol_float_t,Legendre>(poly_indices, abscissas);
  
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
          values[dim][q]
            += normal_vector[dim] * coefficients[hyEdge_dimT * n_shape_fct_ + i] * fct_value;
      }
    }
  }
  
  return values;
}
/*!*********************************************************************************************
 * \brief   Evaluate the function lambda on tensor product points on the boundary
 *
 * \todo not yet implemented
 *
 * \tparam  absc_float_t  Floating type for the abscissa values.
 * \tparam  sizeT         Size of the array of array of abscissas.
 * \param   abscissas     Abscissas of the supporting points.
 * \param   lambda_values The values of the skeletal variable's coefficients.
 * \param   boundary_number number of the boundary on which to evaluate the function.
 **********************************************************************************************/
template<unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg, unsigned int quad_deg, typename lSol_float_t>
template<typename abscissa_float_t, std::size_t sizeT, class input_array_t>
std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>, 1> TimoschenkoBendingBeam<
    hyEdge_dimT,
    space_dim,
    poly_deg,
    quad_deg,
    lSol_float_t>::lambda_values(const std::array<abscissa_float_t, sizeT> &abscissas,
                                   const input_array_t &lambda_values,
                                   const unsigned int boundary_number) const {
  return std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(sizeT)>, 1>();
}
// end of TimoschenkoBendingBeam::bulk_values


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
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3). For this class, all
 * hyperedges are supposed to be uniform (i.e. equal to the unit hypercube). Thus, no geometrical
 * information is needed by this class.
 *
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dim, unsigned int poly_deg,
  unsigned int quad_deg, typename lSol_float_t = double
>
class LengtheningTimoschenkoBendingBeam
{
  public:
    /*!*********************************************************************************************
     * \brief   Return template parameter \c hyEdge_dimT.
     * 
     * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Specify whether advanced functuions are implemented for this class
     * 
     * Advanced functions are numerical_flux_from_rhs, dirichlet_coeffs, and neumann_coeffs.
     * 
     * \todo    Replace this by concept in C++20!
     **********************************************************************************************/    
    static constexpr bool advanced_functions() { return false; }
    /*!*********************************************************************************************
     * \brief   Evaluate amount of global degrees of freedom per hypernode.
     * 
     * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
     * cooperating with this object.
     *
     * \retval  n_dofs        Number of global degrees of freedom per hypernode.
     **********************************************************************************************/
    static constexpr unsigned int n_glob_dofs_per_node()
    { return 2 * space_dim * Hypercube<hyEdge_dimT-1>::pow(poly_deg + 1); }
    
    
    static constexpr unsigned int node_value_dimension() { return space_dim; }
    
    static constexpr unsigned int system_dimension() { return space_dim; }
    
    
  private:
    
    const lSol_float_t tau_;

    const LengtheningBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t> len_beam;
    const TimoschenkoBendingBeam<hyEdge_dimT,space_dim,poly_deg,quad_deg,lSol_float_t> ben_beam;

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
    LengtheningTimoschenkoBendingBeam(const constructor_value_type& tau = 1.)
    : tau_(tau), len_beam(tau), ben_beam(tau)  { }
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
    template <class hyEdgeT>
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT >
    numerical_flux_from_lambda
    (const std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2*hyEdge_dimT >&
      lambda_values, hyEdgeT& hyper_edge ) const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > result, aux;
      
      result = len_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);
      aux = ben_beam.numerical_flux_from_lambda(lambda_values, hyper_edge);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int j = 0; j < n_glob_dofs_per_node(); ++j)
          result[i][j] += aux[i][j];

      return result;
    }

    /*!*********************************************************************************************
     * \brief   Evaluate local contribution to right-hand side vector by global right-hand side.
     *
     * \tparam  GeomT         The geometry type / typename of the considered hyEdge's geometry.
     * \param   geom          The geometry of the considered hyperedge (of typename GeomT).
     * \retval  vec_b         Local part of vector b.
     **********************************************************************************************/
    template < class hyEdgeT >
    std::array< std::array<lSol_float_t, n_glob_dofs_per_node()>, 2*hyEdge_dimT > numerical_flux_from_rhs
    ( const std::array< unsigned int, 2*hyEdge_dimT > & bound_cond, hyEdgeT& hyper_edge )  const
    {
      std::array< std::array<lSol_float_t, n_glob_dofs_per_node()> , 2 * hyEdge_dimT > bdr_values;            
      return bdr_values;
    }
    
    template<typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()>
    bulk_values
    (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
     hyEdgeT& hyper_edge, const lSol_float_t time = 0.) const
    {
      std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>,system_dimension()>
        result, auxiliary;
      result = len_beam.bulk_values(abscissas,lambda_values,hyper_edge);
      auxiliary = ben_beam.bulk_values(abscissas,lambda_values,hyper_edge);

      for (unsigned int i = 0; i < system_dimension(); ++i)
        for (unsigned int j = 0; j < Hypercube<hyEdge_dimT>::pow(sizeT); ++j)
          result[i][j] += auxiliary[i][j];

      return result;
    }

  template < typename abscissa_float_t, std::size_t sizeT, class input_array_t>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1>
  lambda_values
      (const std::array<abscissa_float_t,sizeT>& abscissas, const input_array_t& lambda_values,
       const unsigned int boundary_number) const
  {
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT-1>::pow(sizeT)>,1> result;
    return result;
  }
    
}; // end of class LengtheningTimoschenkoBendingBeam
