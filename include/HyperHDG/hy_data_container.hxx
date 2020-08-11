#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <vector>


/*!*************************************************************************************************
 * \brief   Hypergraph topology based on an input file.
 *
 * \todo    ALL NEW!
 *
 * The topology class File is a set of hyperedges. Each of these tensorial hyperedges is represented
 * by its hypernodes (given within the file). For consistency, it is assumed that the vertices and 
 * the hypernodes are assumed to be given in lexicographical order to ensure that geometry and 
 * topology of all hyperedges fit.
 * 
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is \c hyEdge_index_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < typename data_t,  typename hyEdge_index_t = unsigned int >
class HyDataContainer
{
  private:
    /*!*********************************************************************************************
     * \brief   Domain Info containing all the information of the hypergraph (cf. ReadDomain.hxx).
     **********************************************************************************************/
    std::vector<data_t> data_container;
    
  public:
    /*!*********************************************************************************************
     * \brief   Defines the return value of the class.
     *
     * The \c class \c HyperGraph_Cubic defines the topology of the hypergraph. It "contains" the
     * different hyperedges (that actually are constructed everytime access is needed from e.g. the
     * solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef data_t value_type;
    /*!*********************************************************************************************
     * \brief   Defines the value type of input argument for standard constructor.
     *
     * To receive a very general \c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The \c typedef implements the aspect, that a cubic hypergraph
     * topology is by default constructed by a std::vector that contains amounts of elements in the
     * different dimensions.
     **********************************************************************************************/
    typedef unsigned int constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Construct a topology from a given filename.
     *
     * \param   filename    Name of file containing the information.
     **********************************************************************************************/
    HyDataContainer(const constructor_value_type n_hyEdges)
    : data_container(n_hyEdges) { }
    
    /*!*********************************************************************************************
     * \brief   Get topological hyperedge of given index.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *        -> This is not really a large object, is it? I mean, it only consists of a reference
     *        and an index. I do not really see the advantage of returning a shared pointer. Does
     *        it make any difference, here?
     *
     * This is equivalent to \c get_hyEdge.
     *
     * \param   index           The index of the hyperedge to be returned.
     * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    value_type& operator[](const hyEdge_index_t index)
    { return get_hyEdge(index); }
    /*!*********************************************************************************************
     * \brief   Get topological hyperedge of given index.
     *
     * This is equivalent to \c operator[].
     *
     * \param   index           The index of the hyperedge to be returned.
     * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    value_type& get_hyEdge(const hyEdge_index_t index)
    {
      hy_assert( index < data_container.size() && index >= 0 ,
                 "Index must be non-negative and smaller than " << data_container.size() <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return data_container[index];
    }
}; // end of class File
