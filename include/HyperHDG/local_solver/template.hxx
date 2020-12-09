#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hypercube.hxx>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief Implementation the interface of local solvers needed by enclosing objects in minimal way.
 *
 * The class serves as a template for the implementation of local solvers. It also serves as a
 * minimal input to other classes using local solvers.The local solver uses no degrees of freedom on
 * the nodes and also produces no output.
 **************************************************************************************************/
template <int hyEdge_dimT, typename lSol_float_t>
class Template
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  typedef struct empty_class
  {
  } data_type;
  /*!***********************************************************************************************
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Number of degrees of freedom per hypernode.
   *
   * \note    \c n_dofs_per_nodeT of HyperNodeFactory must be equal to this number.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node() { return 0U; }
  /*!***********************************************************************************************
   * \brief   The dimension of the local system of partial differential equations
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return 0U; }
  /*!***********************************************************************************************
   * \brief   The dimension of the function represented by skeletal unknowns.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 0U; }
  /*!***********************************************************************************************
   * \brief   The local solver as needed by the HDG method
   ************************************************************************************************/
  template <typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& numerical_flux_from_lambda(const SmallMatInT&, SmallMatOutT&) const
  {
    return SmallMatOutT();
  }
  /*!***********************************************************************************************
   * \brief   The values of the local solution in quadrature points of the cell.
   *
   * \retval  An array of array with outer size system_dimension() and inner size number of
   *          quadrature points, containing for each component of the solution the values in each
   *          quadrature point.
   ************************************************************************************************/
  template <typename AbscissaType, std::size_t AbscissaSize, class InputArrayType>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(AbscissaSize)>,
             system_dimension()>
  bulk_values(const std::array<AbscissaType, AbscissaSize>&,
              const InputArrayType&,
              const lSol_float_t = 0.) const
  {
    return std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(AbscissaSize)>,
                      system_dimension()>();
  }
  /*!***********************************************************************************************
   * \brief   The values of the skeletal variable in quadrature points of faces.
   ************************************************************************************************/
  template <typename AbscissaType, std::size_t AbscissaSize, class InputArrayType>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(AbscissaSize)>,
             node_system_dimension()>
  lambda_values(const std::array<AbscissaType, AbscissaSize>&,
                const InputArrayType&,
                unsigned int) const
  {
    return std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT - 1>::pow(AbscissaSize)>,
                      node_system_dimension()>();
  }
};  // end of class Template

}  // end of namespace LocalSolver
