#pragma once // Ensure that file is included only once in a single compilation.

/**
 * \brief Implementation the interface of local solvers needed by enclosing objects in minimal way.
 *
 * The class serves as a template for the implementation of local solvers. It also serves as a
 * minimal input to other classes using local solvers.
 * The local solver uses no degrees of freedom on the nodes and also produces no output.
 *
 * \todo Add functions for naming and structuring VTK output of ode and edge values
 */
template <int hyEdge_dimT, typename lSol_float_t>
class LocalSolverTemplate
{
public:
  /**
   * \brief   Return template parameter \c hyEdge_dimT.
   * 
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
   */
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  
  /**
   * \brief Number of degrees of freedom per hypernode.
   * 
   * \note This number should be equal to \c n_dofs_per_nodeT of HyperNodeFactory.
   */  
  static constexpr unsigned int n_glob_dofs_per_node()
  { return 0U; }

  /**
   * \brief The dimension of the function space for Lagrange multiplies on nodes
   */
  static constexpr unsigned int node_value_dimension()
  { return 0U; }

  /**
   * \brief The dimension of the local system of partial differential equations
   */
  static constexpr unsigned int system_dimension()
  { return 0U; }
  
  /**
   * \brief The local solver as needed by the HDG method
   */
  std::array< std::array<lSol_float_t, 0> , 2 * hyEdge_dimT > numerical_flux_from_lambda
  (const std::array< std::array<lSol_float_t, 0> , 2*hyEdge_dimT >& lambda_values) const
  {
    return std::array< std::array<lSol_float_t, 0> , 2 * hyEdge_dimT > ();
  }

  /**
   * \brief The values of the local solution in quadrature points of the cell
   * 
   * \retval An array of array with outer size system_dimension() and
   * inner size number of quadrature points, containing for each
   * component of the solution the values in each quadrature point.
   * 
   * \todo  How about the geometry of the respective element?
   */
  template<typename AbscissaType, std::size_t AbscissaSize, class InputArrayType>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(AbscissaSize)>,system_dimension()>
  bulk_values (const std::array<AbscissaType,AbscissaSize>& abscissas,
	       const InputArrayType& lambda_values, const lSol_float_t time = 0.) const
    {
      return std::array<std::array<lSol_float_t,
				   Hypercube<hyEdge_dimT>::pow(AbscissaSize)>,
			system_dimension()> ();
    }
};
