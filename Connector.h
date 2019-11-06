/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "TypeDefs.h"
#include <vector>

template <unsigned int connector_dim, unsigned int space_dim>
class Connector_RegularQuad
{
  private:
    std::vector<joint_index_type> joint_indices_;
    std::vector<bool> correct_joint_orientation_;
  public:
    Connector_RegularQuad(const connector_index_type index, const std::vector<unsigned int>& num_elements,
                          const connector_index_type num_of_connectors);
    const std::vector<joint_index_type>& get_joint_indices() const;
    std::vector<double> abs_det_of_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
    std::vector< std::vector<double> > inv_of_transposed_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
};

#endif
