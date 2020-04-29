#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>

#include <array>
#include <vector>

/*!*************************************************************************************************
 * \brief   This class is responsible for the mapping of hypernodes to global degrees of freedom.
 * 
 * This class administrates the mapping of topological information (i.e. a hypernode) to information
 * related to the approximation of solutions to PDEs. Thus, it is capable of mapping an hypernode
 * index to the indices and values of degrees of freedom. It is the only class to manage the access
 * to the global vectors comprising degrees of freedom and is universal to (almost) all kinds of
 * possible equations.
 * 
 * \tparam  n_dofs_per_nodeT  Amount of degrees of freedom associated to an hypernode. This is the
 *                            number of local trial functions for the skeletal variable.
 * \tparam  hyNode_index_t    Unsigned integer type specification. Default is unsigned int.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_dofs_per_nodeT, typename hyNode_index_t = unsigned int >
class HyperNodeFactory
{
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the amount of dofs per node.
     *
     * \retval  n_dofs_per_node     The amount of degrees of freedom per node.
     **********************************************************************************************/
    static constexpr unsigned int n_dofs_per_node() { return n_dofs_per_nodeT; }
  private:
    /*!*********************************************************************************************
     * \brief   Amount of hypernodes within hypergraph.
     * 
     * The number of hypernodes within the considered hypergraph is needed to construct vectors of
     * the correct size, to check whether a vector has the appropriate size, and to check whether a
     * degree of freedom has a valid index.
     **********************************************************************************************/
    const hyNode_index_t n_hyNodes_;
  public:
    /*!*********************************************************************************************
     * 
     * \brief   Construct HyperNodeFactory from total number of hypernodes. 
     * 
     * \param   n_hypernodes   Total number of hypernodes.
     **********************************************************************************************/
    HyperNodeFactory(const hyNode_index_t n_hyNodes) : n_hyNodes_(n_hyNodes) { }
    /*!*********************************************************************************************
     * \brief   Copy constructot for HypernodeFactory.
     * 
     * Function that evaluates the condensed, matrix-free version of the matrix-vector product
     * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method, \f$x\f$ is the
     * vector of parameters to define the skeletal variable \f$\lambda\f$, and \f$y\f$ is the
     * resulting vector, which has the same size as the input vector \f$x\f$. 
     * 
     * \param   hnf                 A \c HyperNodeFactory to be copied.
     **********************************************************************************************/
    HyperNodeFactory(const HyperNodeFactory<n_dofs_per_nodeT>& hnf) : n_hyNodes_(hnf.n_hyNodes_) { }
    /*!*********************************************************************************************
     * \brief   Returns the total amount of hypernodes in the considered hypergraph.
     * 
     * \retval  n_hypernodes        The total amount of hypernodes in the considered hypergraph.
     **********************************************************************************************/
    const hyNode_index_t n_hyNodes() const { return n_hyNodes_; }
    /*!*********************************************************************************************
     * \brief   Returns the total amount of degrees of freedom in the considered hypergraph.
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \retval  n_global_dofs       The total amount of degreees of freedom in the considered
     *                              hypergraph.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int >
    const dof_index_t n_global_dofs() const  { return n_hyNodes_ * n_dofs_per_nodeT; }
    /*!*********************************************************************************************
     * \brief   Calculate global indices of degrees of freedom related to a hypernode.
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \param   hypernode_index     Index of the considered hypernode.
     * \retval  dof_indices         A \c std::array containing the global indices of related degrees
     *                              of freedom.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int >
    std::array<dof_index_t, n_dofs_per_nodeT> get_dof_indices
      (const hyNode_index_t hyNode_index) const
    {
      dof_index_t initial_dof_index = hyNode_index * n_dofs_per_nodeT;
      std::array<dof_index_t, n_dofs_per_nodeT> dof_indices;
      for (unsigned int i = 0; i < n_dofs_per_nodeT; ++i)
        dof_indices[i] = initial_dof_index + i;
      return dof_indices;
    }
    /*!*********************************************************************************************
     * \brief   Calculate index of hypernode holding a special degree of freedom.
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \param   dof_index           Global index of related degree of freedom.
     * \retval  hypernode_index     Index of the considered hypernode.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int >
    hyNode_index_t get_hyNode_from_dof_index ( const dof_index_t dof_index )  const
    { return dof_index / n_dofs_per_nodeT; }
    /*!*********************************************************************************************
     * \brief   Evaluate values of degrees of freedom related to a hypernode.
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \tparam  dof_value_t         Floating point type specification. Default is double.
     * \param   hypernode_index     Index of the considered hypernode.
     * \retval  dof_values          A \c std::array containing the values of related degrees of
     *                              freedom.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int, typename dof_value_t = double >
    std::array<dof_value_t, n_dofs_per_nodeT> get_dof_values
      (const hyNode_index_t hyNode_index, const std::vector<dof_value_t>& global_dof_vector) const
    {
      dof_index_t initial_dof_index = hyNode_index * n_dofs_per_nodeT;
      hy_assert( initial_dof_index >= 0
                   && initial_dof_index + n_dofs_per_nodeT <= global_dof_vector.size() ,
                 "The initial dof index = " << initial_dof_index << ", should be non-negative. " <<
                 " Moreover, the final index = " << initial_dof_index + n_dofs_per_nodeT << " must "
                 << "not exceed the size of the vector of global degrees of freedom." );
      std::array<dof_value_t, n_dofs_per_nodeT> local_dof_values;
      for (unsigned int index = 0; index < n_dofs_per_nodeT; ++index)
        local_dof_values[index] = global_dof_vector[initial_dof_index + index];
      return local_dof_values;
    }
    /*!*********************************************************************************************
     * \brief   Addy different values to values of degrees of freedom related to a hypernode.
     * 
     * Add local values of the \c std::array \c local_dof_vector to the respective values of the
     * \c std::vector \c global_dof_vector comprising all degrees of freedom. 
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \tparam  dof_value_t         Floating point type specification. Default is double.
     * \param   hypernode_index     Index of the considered hypernode.
     * \param   global_dof_vector   \c std::vector containing the values of all degrees of freedom.
     * \param   local_dof_vector    \c std::array containing the local values to be added to the
     *                              global ones.
     * \retval  global_dof_vector   \c std::vector containing the values of all degrees of freedom.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int, typename dof_value_t = double >
    void add_to_dof_values
      (const hyNode_index_t hyNode_index, std::vector<dof_value_t>& global_dof_vector,
       const std::array<dof_value_t, n_dofs_per_nodeT>& local_dof_vector) const
    {
      dof_index_t initial_dof_index = hyNode_index * n_dofs_per_nodeT;
      hy_assert( local_dof_vector.size() == n_dofs_per_nodeT ,
                 "The size of the local dof vector is " << local_dof_vector.size() << ", but should"
                 << " be equal to the amount of local dofs, which is " << n_dofs_per_nodeT << "." );
      hy_assert( initial_dof_index >= 0
                   && initial_dof_index + n_dofs_per_nodeT <= global_dof_vector.size() ,
                 "The initial dof index = " << initial_dof_index << "should be non-negative. " <<
                 "Moreover, the final index = " << initial_dof_index + n_dofs_per_nodeT << " must "
                 << "not exceed the size of the vector of global degrees of freedom." );
      for(unsigned int index = 0; index < n_dofs_per_nodeT; ++index)
        global_dof_vector[initial_dof_index + index] += local_dof_vector[index];
    }
    /*!*********************************************************************************************
     * \brief   Set all values of degrees of freedom of a hypernode to a predefined value.
     * 
     * \tparam  dof_index_t         Unsigned integer type specification. Default is unsigned int.
     * \tparam  dof_value_t         Floating point type specification. Default is double.
     * \param   hypernode_index     Index of the considered hypernode.
     * \param   global_dof_vector   \c std::vector containing the values of all degrees of freedom.
     * \param   value               The future value of related degrees of freedom.
     * \retval  global_dof_vector   \c std::vector containing the values of all degrees of freedom.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int, typename dof_value_t = double >
    void set_dof_values(const hyNode_index_t hyNode_index,
      std::vector<dof_value_t>& global_dof_vector, const dof_value_t value) const
    {
      dof_index_t initial_dof_index = hyNode_index * n_dofs_per_nodeT;
      hy_assert( initial_dof_index >= 0
                   && initial_dof_index + n_dofs_per_nodeT <= global_dof_vector.size() ,
                 "The initial dof index = " << initial_dof_index << "should be non-negative. " <<
                 "Moreover, the final index = " << initial_dof_index + n_dofs_per_nodeT << " must "
                 << "not exceed the size of the vector of global degrees of freedom." );
      for(unsigned int index = 0; index < n_dofs_per_nodeT; ++index)
        global_dof_vector[initial_dof_index + index] = value;
    }
}; // end of class HyperNodeFactory
