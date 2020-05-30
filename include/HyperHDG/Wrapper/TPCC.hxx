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

/*!*************************************************************************************************
 * \todo  Remove this, when submodule works without it!
 **************************************************************************************************/
#include <ostream>               // TO BE REMOVED!

#include <tpcc/lexicographic.h>  // Submodule which is wrapped by this file!
#include <HyperHDG/HyAssert.hxx>
#include <array>


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

/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t = unsigned int >
tpcc_t < hyEdge_dim, space_dim, index_t >
create_tpcc( const std::array<unsigned int,space_dim>& dimensions )
{ 
  static_assert( space_dim >= hyEdge_dim , "Hypercube dim must not be bigger than spatial dim!");
  return TPCC::Lexicographic<space_dim,hyEdge_dim,index_t>(dimensions);
}
/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
index_t n_elements( const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc ) { return tpcc.size(); }
/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
tpcc_elem_t < hyEdge_dim, space_dim >
get_element( const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc, const index_t index )
{
  hy_assert( index < n_elements(tpcc) , "Index " << index << " must not be bigger than the TPCC "
             << "size, which is " << n_elements(tpcc) << "." );
  return tpcc.operator[](index);
}
/*!*************************************************************************************************
 * \brief   Return index of given element within TPCC.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t >
index_t get_element_index
(  const tpcc_t<hyEdge_dim,space_dim,index_t>& tpcc, const tpcc_elem_t<hyEdge_dim,space_dim>& elem )
{ 
  index_t index = tpcc.index(elem);
  hy_assert( index < n_elements(tpcc) , "Returned index is larger than number of elements!" );
  return index;
}
/*!*************************************************************************************************
 * \brief   Return i-th element facet.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim >
tpcc_elem_t < hyEdge_dim-1, space_dim >
get_element_facet ( const tpcc_elem_t<hyEdge_dim, space_dim>& element, const unsigned int index )
{ 
  hy_assert ( index < 2 * hyEdge_dim , "An element hast only " << 2 * hyEdge_dim << "facets. "
              << "You requested number " << index << "." );
  return element.facet(index);
}