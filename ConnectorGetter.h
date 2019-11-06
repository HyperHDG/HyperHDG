/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef CONNECTOR_GETTER_H
#define CONNECTOR_GETTER_H

#include "TypeDefs.h"
#include "Connector.h"
#include <vector>

/*
template <class AbstractConnector>
class ConnectorGetter
{
  private:
    vector<AbstractConnector> connectors;
  public:
    AbstractConnector& get_connector(const unsigned int index);
    unsigned int num_of_connectors();
};
*/
template <unsigned int connector_dim, unsigned int space_dim>
class ConnectorGetter_RegularQuad
{
  private:
    std::vector<unsigned int> num_elements_;
    connector_index_type num_of_connectors_;
  public:
    ConnectorGetter_RegularQuad(const unsigned int num_of_elem_in_x_dir,
                                const unsigned int num_of_elem_in_y_dir = 0,
                                const unsigned int num_of_elem_in_z_dir = 0);
    ConnectorGetter_RegularQuad
      (const ConnectorGetter_RegularQuad<connector_dim,space_dim>& other);
    
    const Connector_RegularQuad<connector_dim,space_dim> 
      get_connector(const connector_index_type index) const;
    const connector_index_type num_of_connectors() const;
};

#endif
