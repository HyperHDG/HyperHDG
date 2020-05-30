/*!*************************************************************************************************
 * \file    TPP.hxx
 * \brief   This file provides the functions of the submodule TPP.
 *
 * \todo    The external module does not provide functions that return the template parameters. Why?
 *          If these functins will be provided, we should use them to do consistency checks!
 *
 * This is a wrapper file to provide TPP (tensor product chain complex) based functions. Since this
 * is an external submodule utilized by hyperHDG, this file provides appropriate wrapper functions.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/

#pragma once // Ensure that file is included only once in a single compilation.

#include <submodules/tensor_product_chain_complex.git/include/tpcc/lexicographic.h>
#include <HyperHDG/HyAssert.hxx>
#include <array>


/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex.
 **************************************************************************************************/
template < unsigned int hyEdge_dim, unsigned int space_dim, typename index_t = unsigned int >
auto create_tpcc(const std::array<unsigned int, space_dim>& dimensions)
{ 
  static_assert( space_dim >= hyEdge_dim , "Hypercube dim must not be bigger than spatial dim!");
  return Lexicographic<space_dim,hyEdge_dim,index_t>(dimensions);
}
/*!*************************************************************************************************
 * \brief   Return the element of given index the TPP.
 **************************************************************************************************/
auto get_element(const auto& tpcc, const hyEdge_index_t index)
{
  hy_assert( index < tpcc.size() ,
             "Index " + index + " must not be bigger than the TPP size " + tpp.size() + "." );
  return tpcc.operator[](index);
}
/*!*************************************************************************************************
 * \brief   Return index of given element within TPP.
 **************************************************************************************************/
auto get_element_index(const auto& tpcc, const auto& element)
{ return tpcc.index(element); }
/*!*************************************************************************************************
 * \brief   Return i-th element facet.
 **************************************************************************************************/
auto get_element_facet(const auto& element, const unsigned int index)
{ return hyEdge.facet(index); }