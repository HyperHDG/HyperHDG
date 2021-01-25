/*!*************************************************************************************************
 * \file    tpcc.hxx
 * \brief   This file provides the functions of the submodule TPCC.
 *
 * This is a wrapper file to provide TPCC (tensor product chain complex) based functions. Since this
 * is an external submodule utilized by hyperHDG, this file provides appropriate wrapper functions.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/

#pragma once  // Ensure that file is included only once in a single compilation.

#include <tpcc/lexicographic.h>  // Submodule which is wrapped by this file!
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>

namespace Wrapper
{
// -------------------------------------------------------------------------------------------------
// Wrapper classes for mathematical helper functions.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Binomial coefficient of unsigned integers.
 **************************************************************************************************/
unsigned int binomial(unsigned int n, unsigned int k)
{
  return TPCC::binomial<unsigned int>(n, k);
}

// -------------------------------------------------------------------------------------------------
// Wrapper classes for tensor product chain complex and its elements.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Type of a tensor product chain complex.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          typename TPCC::boundaries bndT = TPCC::boundaries::both,
          typename index_t = unsigned int>
using tpcc_t =
  TPCC::Lexicographic<space_dim, hyEdge_dim, bndT, index_t, unsigned int, unsigned int>;
/*!*************************************************************************************************
 * \brief   Type of an element of a tensor product chain complex.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim>
using tpcc_elem_t = TPCC::Element<space_dim, hyEdge_dim, unsigned int, unsigned int>;

// -------------------------------------------------------------------------------------------------
// Functions considering full tensor product chain complexes.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          typename TPCC::boundaries bndT = TPCC::boundaries::both,
          typename index_t = unsigned int>
tpcc_t<hyEdge_dim, space_dim, bndT, index_t> create_tpcc(const SmallVec<space_dim, index_t>& vec)
{
  static_assert(space_dim >= hyEdge_dim, "Hypercube dim must not be bigger than spatial dim!");
  return TPCC::Lexicographic<space_dim, hyEdge_dim, bndT, index_t, unsigned int, unsigned int>(
    vec.data());
}
/*!*************************************************************************************************
 * \brief   Create a tensor product chain complex associated to the facets.
 **************************************************************************************************/
template <typename TPCC::boundaries bndT>
auto tpcc_faces(const auto& elements)
{
  return elements.template boundary<bndT>();
}
/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
template <typename index_t = unsigned int>
index_t n_elements(const auto& tpcc)
{
  return tpcc.size();
}

// -------------------------------------------------------------------------------------------------
// Functions related to single elements.
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Return the element of given index the TPCC.
 **************************************************************************************************/
auto get_element(const auto& tpcc, const auto index)
{
  hy_assert(index < tpcc.size(), "Index " << index << " must not be bigger than the TPCC "
                                          << "size, which is " << tpcc.size() << ".");
  return tpcc.operator[](index);
}
/*!*************************************************************************************************
 * \brief   Return index of given element within TPCC.
 **************************************************************************************************/
template <typename index_t = unsigned int>
index_t get_index(const auto& tpcc, const auto& elem)
{
  index_t index = tpcc.index(elem);
  hy_assert(index < tpcc.size(), "Returned index is larger than number of elements!");
  return index;
}
/*!*************************************************************************************************
 * \brief   Return i-th element facet.
 **************************************************************************************************/
auto get_face(const auto& elem, const unsigned int index)
{
  hy_assert(index < elem.n_facets(), "An element hast only " << elem.n_facets() << "facets. "
                                                             << "You requested number " << index
                                                             << ".");
  return elem.facet(index);
}
/*!*************************************************************************************************
 * \brief   Determine the orientation of an element.
 **************************************************************************************************/
unsigned int elem_orientation(const auto& elem)
{
  return elem.direction_index();
}
/*!*************************************************************************************************
 * \brief   Return the index-th orthonormal direction to element.
 **************************************************************************************************/
unsigned int exterior_direction(const auto& elem, const unsigned int index)
{
  hy_assert(index < elem.n_val - elem.k_val,
            "There are only " << elem.n_val - elem.k_val << " exterior directions.");
  unsigned int acr_dir = elem.across_direction(index);
  hy_assert(acr_dir < elem.n_val,
            "Exterior direction must be smaller than amount of spatial dimensions!");
  return acr_dir;
}
/*!*************************************************************************************************
 * \brief   Return coordinate value with respect to index-th orthonormal direction of element.
 **************************************************************************************************/
unsigned int exterior_coordinate(const auto& elem, const unsigned int index)
{
  hy_assert(index < elem.n_val - elem.k_val,
            "There are only " << elem.n_val - elem.k_val << " exterior directions.");
  return elem.across_coordinate(index);
}

}  // end of namespace Wrapper
