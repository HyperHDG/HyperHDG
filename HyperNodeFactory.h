/*!*************************************************************************************************
 * @brief   This class is responsible for the mapping of hypernodes to global degrees of freedom.
 * 
 * This class administrates the mapping of topological information (i.e. a hypernode) to information
 * related to the approximation of solutions to PDEs. Thus, it is capable of mapping an hypernode
 * index to the indices and values of degrees of freedom. It is the only class to manage the access
 * to the global vectors comprising degrees of freedom and is universal to (almost) all kinds of
 * possible equations.
 * 
 * @tparam  n_dofs_per_node   Amount of degrees of freedom associated to an hypernode. This is the
 *                            number of local trial functions for the skeletal variable.
 * 
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/

#ifndef HYPERNODEFACTORY_H
#define HYPERNODEFACTORY_H

#include "TypeDefs.h"
#include <array>
#include <vector>

template <unsigned int n_dofs_per_node>
class HyperNodeFactory
{
  private:
    /*!*********************************************************************************************
     * @brief   Amount of hypernodes within hypergraph.
     * 
     * The number of hypernodes within the considered hypergraph is needed to construct vectors of
     * the correct size, to check whether a vector has the appropriate size, and to check whether a
     * degree of freedom has a valid index.
     **********************************************************************************************/
    const hypernode_index_type num_of_hypernodes_;
  public:
    /*!*********************************************************************************************
     * 
     * @brief   Construct HyperNodeFactory from total number of hypernodes. 
     * 
     * @param   num_of_hypernodes   Total number of hypernodes.
     **********************************************************************************************/
    HyperNodeFactory(const hypernode_index_type num_of_hypernodes);
    /*!*********************************************************************************************
     * @brief   Copy constructot for HypernodeFactory.
     * 
     * Function that evaluates the condensed, matrix-free version of the matrix-vector product
     * @f$A x = y@f$, where @f$A@f$ is the condensed matrix of the LDG-H method, @f$x@f$ is the
     * vector of parameters to define the skeletal variable @f$\lambda@f$, and @f$y@f$ is the
     * resulting vector, which has the same size as the input vector @f$x@f$. 
     * 
     * @param   other               A @c HyperNodeFactory to be copied.
     **********************************************************************************************/
    HyperNodeFactory(const HyperNodeFactory<n_dofs_per_node>& other);
    /*!*********************************************************************************************
     * @brief   Returns the total amount of hypernodes in the considered hypergraph.
     * 
     * @retval  num_of_hypernodes   The total amount of hypernodes in the considered hypergraph.
     **********************************************************************************************/
    const hypernode_index_type num_of_hypernodes() const;
    /*!*********************************************************************************************
     * @brief   Returns the total amount of degrees of freedom in the considered hypergraph.
     * 
     * @retval  num_of_global_dofs  The total amount of degreees of freedom in the considered
     *                              hypergraph.
     **********************************************************************************************/
    const dof_index_type num_of_global_dofs() const;
    /*!*********************************************************************************************
     * @brief   Calculate global indices of degrees of freedom related to a hypernode.
     * 
     * @param   hypernode_index     Index of the considered hypernode.
     * @retval  dof_indices         A @c std::array containing the global indices of related degrees
     *                              of freedom.
     **********************************************************************************************/
    std::array<dof_index_type, n_dofs_per_node> get_dof_indices
      (const hypernode_index_type hypernode_index) const;
    /*!*********************************************************************************************
     * @brief   Evaluate values of degrees of freedom related to a hypernode.
     * 
     * @param   hypernode_index     Index of the considered hypernode.
     * @retval  dof_values          A @c std::array containing the values of related degrees of
     *                              freedom.
     **********************************************************************************************/
    std::array<dof_value_type, n_dofs_per_node> get_dof_values
      (const hypernode_index_type hypernode_index,
       const std::vector<dof_value_type>& global_dof_vector) const;
    /*!*********************************************************************************************
     * @brief   Addy different values to values of degrees of freedom related to a hypernode.
     * 
     * Add local values of the @c std::array @c local_dof_vector to the respective values of the
     * @c std::vector @c global_dof_vector comprising all degrees of freedom. 
     * 
     * @param   hypernode_index     Index of the considered hypernode.
     * @param   global_dof_vector   @c std::vector containing the values of all degrees of freedom.
     * @param   local_dof_vector    @c std::array containing the local values to be added to the
     *                              global ones.
     * @retval  global_dof_vector   @c std::vector containing the values of all degrees of freedom.
     **********************************************************************************************/
    void add_to_dof_values
      (const hypernode_index_type hypernode_index, std::vector<dof_value_type>& global_dof_vector,
       const std::array<dof_value_type, n_dofs_per_node>& local_dof_vector) const;
    /*!*********************************************************************************************
     * @brief   Set all values of degrees of freedom of a hypernode to a predefined value.
     * 
     * @param   hypernode_index     Index of the considered hypernode.
     * @param   global_dof_vector   @c std::vector containing the values of all degrees of freedom.
     * @param   value               The future value of related degrees of freedom.
     * @retval  global_dof_vector   @c std::vector containing the values of all degrees of freedom.
     **********************************************************************************************/
    void set_dof_values(const hypernode_index_type hypernode_index,
      std::vector<dof_value_type>& global_dof_vector, const dof_value_type value) const;
    
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the amount of dofs per node.
     *
     * @retval  n_dofs_per_node     The amount of degrees of freedom per node.
     **********************************************************************************************/
    static constexpr unsigned int n_dof_per_node() { return n_dofs_per_node; };
};

/*!*************************************************************************************************
 * @brief   Calculate the amount of local degrees of freedom of a hypernode at compile time.
 * 
 * Naive implementation without math packages of
 * "amount = solution_dim * (poly_degree ^ (hyperedge_dim - 1))"!
 * 
 * \todo    Find better location, since this can vary with different test functions (P_k, Q_k, ...).
 * 
 * @param   hyperedge_dim       The dimension of a hyperedge (1 for graphs).
 * @param   poly_degree         The local polynomial degree of test functions.
 * @param   solution_dim        The dimension of the solution (1 for scalar equations).
 * @retval  n_dofs_per_node     The amount of degrees of freedom per hypernode.
 * 
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/
constexpr const unsigned int compute_n_dofs_per_node(const unsigned int hyperedge_dim,
  const unsigned int poly_degree, const unsigned int solution_dim = 1)
{
  unsigned int amount = 1;
  for (unsigned int iteration = 0; iteration < hyperedge_dim - 1; ++ iteration)
    amount *= poly_degree + 1;
  amount *= solution_dim;
  return amount;
}

#endif
