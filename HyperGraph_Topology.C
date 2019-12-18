/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperGraph_Topology.h"
#include <cassert>

using namespace std;
using namespace Topology;
#include "HyperGraph_Topology.inst"


template <unsigned int hyperedge_dim, unsigned int space_dim>
HyperGraph_Cubic<hyperedge_dim,space_dim>::
HyperGraph_Cubic(const array<unsigned int, space_dim>& num_elements)
: num_elements_(num_elements)
{
  static_assert( hyperedge_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= hyperedge_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set num_of_hyperedges_
  num_of_hyperedges_ = 1;
  if ( hyperedge_dim == space_dim )
    for (unsigned int dim = 0; dim < space_dim; ++dim)  num_of_hyperedges_ *= num_elements[dim];
  else if ( hyperedge_dim == space_dim - 1 )
  {
    num_of_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      num_of_hyperedges_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 2 )
  {
    num_of_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      num_of_hyperedges_ += helper;
    }
  }
  else  assert( 0 == 1 );  
  assert( num_of_hyperedges_ > 0 );
  
  // Set num_of_hypernodes
  num_of_hypernodes_ = 1;
  if (hyperedge_dim == 1)
  {
    num_of_hypernodes_ *= num_elements[0] + 1;
    if (space_dim > 1)  num_of_hypernodes_ *= num_elements[1] + 1;
    if (space_dim > 2)  num_of_hypernodes_ *= num_elements[2] + 1;
  }
  else if ( hyperedge_dim == space_dim )
  {
    num_of_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      num_of_hypernodes_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 1 )
  {
    num_of_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      num_of_hypernodes_ += helper;
    }
  }
  else  assert( 0 == 1 );
  assert( num_of_hypernodes_ > 0 );
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
HyperGraph_Cubic<hyperedge_dim,space_dim>::
HyperGraph_Cubic(const std::vector<int>& num_elements)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim) num_elements_[dim] = num_elements[dim];
  
  static_assert( hyperedge_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= hyperedge_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set num_of_hyperedges_
  num_of_hyperedges_ = 1;
  if ( hyperedge_dim == space_dim )
    for (unsigned int dim = 0; dim < space_dim; ++dim)  num_of_hyperedges_ *= num_elements[dim];
  else if ( hyperedge_dim == space_dim - 1 )
  {
    num_of_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      num_of_hyperedges_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 2 )
  {
    num_of_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      num_of_hyperedges_ += helper;
    }
  }
  else  assert( 0 == 1 );  
  assert( num_of_hyperedges_ > 0 );
  
  // Set num_of_hypernodes
  num_of_hypernodes_ = 1;
  if (hyperedge_dim == 1)
  {
    num_of_hypernodes_ *= num_elements[0] + 1;
    if (space_dim > 1)  num_of_hypernodes_ *= num_elements[1] + 1;
    if (space_dim > 2)  num_of_hypernodes_ *= num_elements[2] + 1;
  }
  else if ( hyperedge_dim == space_dim )
  {
    num_of_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      num_of_hypernodes_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 1 )
  {
    num_of_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      num_of_hypernodes_ += helper;
    }
  }
  else  assert( 0 == 1 );
  assert( num_of_hypernodes_ > 0 );
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
HyperGraph_Cubic<hyperedge_dim,space_dim>::
HyperGraph_Cubic(const HyperGraph_Cubic<hyperedge_dim,space_dim>& other)
: num_elements_(other.num_elements_), num_of_hyperedges_(other.num_of_hyperedges_),
  num_of_hypernodes_(other.num_of_hypernodes_) { }


template <unsigned int hyperedge_dim, unsigned int space_dim>
const HyperEdge_Cubic<hyperedge_dim, space_dim>
HyperGraph_Cubic<hyperedge_dim,space_dim>::
get_hyperedge(const hyperedge_index_type index) const
{
  assert ( index < num_of_hyperedges_ );
  return HyperEdge_Cubic<hyperedge_dim,space_dim>(index, num_elements_);
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const array<unsigned int, space_dim>& 
HyperGraph_Cubic<hyperedge_dim,space_dim>::
num_elements() const
{
  return num_elements_;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const hyperedge_index_type
HyperGraph_Cubic<hyperedge_dim,space_dim>::
num_of_hyperedges() const
{
  return num_of_hyperedges_;
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
const hypernode_index_type
HyperGraph_Cubic<hyperedge_dim,space_dim>::
num_of_hypernodes() const
{
  return num_of_hypernodes_;
}
