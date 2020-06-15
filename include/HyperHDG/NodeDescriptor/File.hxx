#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/ReadDomain.hxx>

#include <array>
#include <string>

/*!*************************************************************************************************
 * \brief   A namespace containing different classes describing local hypernode properties.
 *
 * \todo    All doxygen!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
namespace NodeDescriptor
{  

/*!*************************************************************************************************
 * \brief   Hypergraph topology based on an input file.
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
template
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyEdge_index_t = unsigned int,
  typename hyNode_index_t = hyEdge_index_t
>
class File
{
  
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
    public:
      static constexpr unsigned int n_hyNodes() { return 2 * hyEdge_dimT; }
      /*!*******************************************************************************************
       * \brief   Returns dimension of the hyperedge.
       ********************************************************************************************/
      static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
      /*!*******************************************************************************************
       * \brief   Returns dimension of the surrounding space.
       ********************************************************************************************/
      static constexpr unsigned int space_dim() { return space_dimT; }
    private:
      /*!*******************************************************************************************
       * \brief   Reference to parent hypergraph.
       ********************************************************************************************/
      const File& hyGraph_topology_;
      /*!*******************************************************************************************
       * \brief   Index of the hyperedge within the hypergraph
       ********************************************************************************************/
      const hyEdge_index_t index_;
    public:
      /*!*******************************************************************************************
       * \brief   Construct hyperedge from hypergraph and index.
       ********************************************************************************************/
      hyEdge ( const File& hyGraph_topology, const hyEdge_index_t index )
      : hyGraph_topology_(hyGraph_topology), index_(index) { }
      /*!*******************************************************************************************
       * \brief   Return hypernodes of a hyperedge.
       ********************************************************************************************/
      const std::array<hyNode_index_t, n_hyNodes()>& get_hyFaces_types() const
      { return hyGraph_topology_.domain_info_.hyFaces_hyEdge[index_]; }
      const unsigned int operator[](const unsigned int index) const 
      { return hyGraph_topology_.domain_info_.hyFaces_hyEdge[index_][index]; }
  }; // end of class hyEdge
  
  public:
    /*!*********************************************************************************************
     * \brief   Return local dimension of the hypergraph's hyperedges.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Return dimension of the surrounding space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; }
  
  private:
    /*!*********************************************************************************************
     * \brief   Domain Info containing all the information of the hypergraph (cf. ReadDomain.hxx).
     **********************************************************************************************/
    const DomainInfo<hyEdge_dimT,space_dimT>& domain_info_;
    
  public:
    /*!*********************************************************************************************
     * \brief   Defines the return value of the class.
     *
     * The \c class \c HyperGraph_Cubic defines the topology of the hypergraph. It "contains" the
     * different hyperedges (that actually are constructed everytime access is needed from e.g. the
     * solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef hyEdge value_type;
    /*!*********************************************************************************************
     * \brief   Defines the value type of input argument for standard constructor.
     *
     * \todo    Doxygen
     **********************************************************************************************/
    typedef Topology::File<hyEdge_dimT,space_dimT> constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    File(const constructor_value_type& topology) : domain_info_(topology.domain_info()) { }
    /*!*********************************************************************************************
     * \brief   Copy constructor.
     **********************************************************************************************/
    File(const File<hyEdge_dimT,space_dimT>& other)
    : domain_info_(other.domain_info) { }
    
    /*!*********************************************************************************************
     * \brief   Get geometrical hyperedge of given index.
     *
     * \todo    Doxygen
     *
     * \param   index       The index of the hyperedge to be returned.
     * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    value_type operator[](const hyEdge_index_t index) const
    { return get_hyEdge(index); }
    /*!*********************************************************************************************
     * \brief   Get geometrical hyperedge of given index.
     *
     * \todo    Doxygen
     *
     * \param   index       The index of the hyperedge to be returned.
     * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    value_type get_hyEdge(const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
}; // end of class File

} // end of namespace Topology
