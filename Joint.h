/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef JOINT_H
#define JOINT_H

#include "TypeDefs.h"
#include <vector>

class Joint_RegularQuad
{
  private:
    const unsigned int initial_dof_index_, amount_of_local_dofs_; // Im Konstruktor uebergeben und daher bekannt
  public:
    Joint_RegularQuad(const unsigned int index, const unsigned int amount_of_local_dofs);
    std::vector<dof_index_type> get_dof_indices() const;
    std::vector<dof_value_type> get_dof_values(const std::vector<dof_value_type>& global_dof_vector) const;
    void add_to_dof_values(std::vector<dof_value_type>& global_dof_vector, const std::vector<dof_value_type>& local_dof_vector) const;
    void set_dof_values(std::vector<dof_value_type>& global_dof_vector, const double value) const;
    std::vector<double> abs_det_of_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
};

#endif
