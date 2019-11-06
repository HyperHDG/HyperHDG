/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "ConnectorGetter.h"
#include <cassert>

using namespace std;


template class ConnectorGetter_RegularQuad<1,1>;
template class ConnectorGetter_RegularQuad<1,2>;
template class ConnectorGetter_RegularQuad<1,3>;
template class ConnectorGetter_RegularQuad<2,2>;
template class ConnectorGetter_RegularQuad<2,3>;
template class ConnectorGetter_RegularQuad<3,3>;


template <unsigned int connector_dim, unsigned int space_dim>
ConnectorGetter_RegularQuad<connector_dim,space_dim>::
ConnectorGetter_RegularQuad(const unsigned int num_of_elem_in_x_dir,
                            const unsigned int num_of_elem_in_y_dir,
                            const unsigned int num_of_elem_in_z_dir)
{
  static_assert( connector_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= connector_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
  
  num_elements_.resize(space_dim);
  num_elements_[0] = num_of_elem_in_x_dir;
  if (space_dim > 1)  num_elements_[1] = num_of_elem_in_y_dir;
  if (space_dim > 2)  num_elements_[2] = num_of_elem_in_z_dir;
    
  num_of_connectors_ = 1;
  if ( connector_dim == space_dim )
    for (unsigned int dim = 0; dim < space_dim; ++dim)  num_of_connectors_ *= num_elements_[dim];
  else if ( connector_dim == space_dim - 1 )
  {
    num_of_connectors_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements_[dim_n] + 1;
        else                 helper *= num_elements_[dim_n];
      num_of_connectors_ += helper;
    }
  }
  else if ( connector_dim == space_dim - 2 )
  {
    num_of_connectors_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements_[dim_n];
        else                 helper *= num_elements_[dim_n] + 1;
      num_of_connectors_ += helper;
    }
  }
  else  assert( 0 == 1 );
  
  assert( num_of_connectors_ > 0 );
}


template <unsigned int connector_dim, unsigned int space_dim>
ConnectorGetter_RegularQuad<connector_dim,space_dim>::
ConnectorGetter_RegularQuad(const ConnectorGetter_RegularQuad<connector_dim,space_dim>& other)
: num_elements_(other.num_elements_), num_of_connectors_(other.num_of_connectors_) { }


template <unsigned int connector_dim, unsigned int space_dim>
const Connector_RegularQuad<connector_dim,space_dim>
ConnectorGetter_RegularQuad<connector_dim,space_dim>::
get_connector(const connector_index_type index) const
{
  assert ( index < num_of_connectors_ );
  return Connector_RegularQuad<connector_dim,space_dim>(index, num_elements_, num_of_connectors_);
}


template <unsigned int connector_dim, unsigned int space_dim>
const connector_index_type
ConnectorGetter_RegularQuad<connector_dim,space_dim>::
num_of_connectors() const
{
  return num_of_connectors_;
}
