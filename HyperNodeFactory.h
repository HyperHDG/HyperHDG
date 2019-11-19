/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef HYPERNODEFACTORY_H
#define HYPERNODEFACTORY_H

#include "TypeDefs.h"
#include <array>
#include <vector>

// Naive implementation without math packages of "amount = (polynomial_degree) ^ (hyperedge_dim - 1)"!
constexpr const unsigned int local_dof_amount_node(const unsigned int hyperedge_dim, const unsigned int polynomial_degree)
{
  unsigned int amount = 1;
  for (unsigned int iteration = 0; iteration < hyperedge_dim - 1; ++ iteration)  amount *= polynomial_degree + 1;
  return amount;
}

template <unsigned int amount_of_local_dofs>
class HyperNodeFactory
{
  private:
    const hypernode_index_type num_of_hypernodes_;
  public:
    HyperNodeFactory(const hypernode_index_type num_of_hypernodes);
    HyperNodeFactory(const HyperNodeFactory<amount_of_local_dofs>& other);
    
    const hypernode_index_type num_of_hypernodes() const;
    const dof_index_type num_of_global_dofs() const;
    
    std::array<dof_index_type, amount_of_local_dofs> get_dof_indices(const hypernode_index_type hypernode_index) const;
    std::array<dof_value_type, amount_of_local_dofs> get_dof_values(const hypernode_index_type hypernode_index, const std::vector<dof_value_type>& global_dof_vector) const;
    void add_to_dof_values(const hypernode_index_type hypernode_index, std::vector<dof_value_type>& global_dof_vector, const std::array<dof_value_type, amount_of_local_dofs>& local_dof_vector) const;
    void set_dof_values(const hypernode_index_type hypernode_index, std::vector<dof_value_type>& global_dof_vector, const dof_value_type value) const;
};

#endif
