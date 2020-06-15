/*!*************************************************************************************************
 * \file    TPCC.hxx
 * \brief   This file provides the functions of the submodule TPCC.
 *
 * \todo    The external module does not provide functions that return the template parameters. Why?
 *          If these functins will be provided, we should use them to do consistency checks!
 *
 * This is a wrapper file to provide TPCC (tensor product chain complex) based functions. Since this
 * is an external submodule utilized by hyperHDG, this file provides appropriate wrapper functions.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/

#pragma once // Ensure that file is included only once in a single compilation.

#include <tpcc/lexicographic.h>  // Submodule which is wrapped by this file!
#include <HyperHDG/HyAssert.hxx>
#include <array>

// -------------------------------------------------------------------------------------------------
// Wrapper classes for mathematical helper functions.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Binomial coefficient of unsigned integers.
 **************************************************************************************************/
unsigned int binomial(unsigned int n, unsigned int k)  { return TPCC::binomial<unsigned int>(n,k); }

// -------------------------------------------------------------------------------------------------
// Wrapper classes for tensor product chain complex and its elements.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Type of a tensor product chain complex.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t = unsigned int >
using tpcc_t = TPCC::Lexicographic<space_dim, hyEdge_dim, index_t, unsigned int, unsigned int>;
/*!*************************************************************************************************
 * \brief   Type of an element of a tensor product chain complex.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim >
using tpcc_elem_t = TPCC::Element<space_dim, hyEdge_dim, unsigned int, unsigned int>;

// -------------------------------------------------------------------------------------------------
// Functions considering full tensor product chain complexes.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t = unsigned int >
tpcc_t < hyEdge_dim, space_dim, index_t >
create_tpcc( const std::array<index_t, space_dim>& dimensions )
{ 
  static_assert( space_dim >= hyEdge_dim , "Hypercube dim must not be bigger than spatial dim!");
  return TPCC::Lexicographic<space_dim,hyEdge_dim,index_t,unsigned int, unsigned int>(dimensions);
}
/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex associated to the facets.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
tpcc_t < hyEdge_dim-1, space_dim, index_t >
tpcc_faces ( const tpcc_t < hyEdge_dim, space_dim, index_t >& elements )
{ return elements.boundary(); }
/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
index_t n_elements( const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc ) { return tpcc.size(); }

// -------------------------------------------------------------------------------------------------
// Functions related to single elements.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
tpcc_elem_t < hyEdge_dim, space_dim >
get_element( const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc, const index_t index )
{
  index_t n_elements_ = n_elements<hyEdge_dim,space_dim,index_t>(tpcc);
  hy_assert( index < n_elements_ , "Index " << index << " must not be bigger than the TPCC "
             << "size, which is " << n_elements_ << "." );
  return tpcc.operator[](index);
}
/*!*************************************************************************************************
 * \brief   Return index of given element within TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
index_t get_index
(  const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc, const tpcc_elem_t<hyEdge_dim,space_dim>& elem )
{ 
  index_t index = tpcc.index(elem);
  index_t n_elements_ = n_elements<hyEdge_dim,space_dim,index_t>(tpcc);
  hy_assert( index < n_elements_ , "Returned index is larger than number of elements!" );
  return index;
}
/*!*************************************************************************************************
 * \brief   Return i-th element facet.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim >
tpcc_elem_t < hyEdge_dim-1, space_dim >
get_face ( const tpcc_elem_t<hyEdge_dim, space_dim>& element, const unsigned int index )
{ 
  hy_assert ( index < 2 * hyEdge_dim , "An element hast only " << 2 * hyEdge_dim << "facets. "
              << "You requested number " << index << "." );
  return element.facet(index);
}

template < unsigned int hyEdge_dim, unsigned int space_dim >
unsigned int elem_orientation ( const tpcc_elem_t<hyEdge_dim,space_dim>& elem )
{ return elem.direction_index(); }


template < unsigned int hyEdge_dim, unsigned int space_dim >
unsigned int exterior_direction
( const tpcc_elem_t<hyEdge_dim,space_dim>& elem, const unsigned int index )
{ 
  hy_assert( index < space_dim - hyEdge_dim ,
             "There are only " << space_dim - hyEdge_dim << " exterior directions." );
  unsigned int acr_dir = elem.across_direction(index);
  hy_assert( acr_dir < space_dim ,
             "Exterior direction must be smaller than amount of spatial dimensions!" );
  return acr_dir;
}


template < unsigned int hyEdge_dim, unsigned int space_dim >
unsigned int exterior_coordinate
( const tpcc_elem_t<hyEdge_dim,space_dim>& elem, const unsigned int index )
{ 
  hy_assert( index < space_dim - hyEdge_dim ,
             "There are only " << space_dim - hyEdge_dim << " exterior directions." );
  return elem.across_coordinate(index);
}
