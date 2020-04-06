#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/DenseLA.hxx>
#include <tensor_mapping.hxx>

#include <HyperHDG/Mapping/Linear.hxx>

#include <array>

/*!*************************************************************************************************
 * \brief   A namespace containing classes describing hypergraph geometries.
 *
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain \f$\Omega\f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
namespace Geometry
{

/*!*************************************************************************************************
 * \brief   Hypergraph geometry based on an input file.
 *
 * \todo    Some functions of this class only work for one-dimensional hyperedges. These need to be
 *          generalized to arbitrary domensions of hyperedges.
 *
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template 
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, typename pt_coord_t = double,
  typename mapping_t = Parallelopipedon<hyEdge_dimT,space_dimT,pt_coord_t>,
  typename hyEdge_index_t = unsigned int
>
class File
{
  
  /*!***********************************************************************************************
   * \brief   Definition of the geometry of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
    public:
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
      const File& hyGraph_geometry_;
      /*!*******************************************************************************************
       * \brief   Index of the hyperedge within the hypergraph
       ********************************************************************************************/
      const hyEdge_index_t index_;
      /*!*******************************************************************************************
       * \brief   Hold an instance of a mapping type to be able to calculate normals and so on.
       ********************************************************************************************/
//      mapping_t mapping;
    public:
      /*!*******************************************************************************************
       * \brief   Construct hyperedge from hypergraph and index.
       ********************************************************************************************/
      hyEdge ( const File& hyGraph_geometry, const hyEdge_index_t index )
      : hyGraph_geometry_(hyGraph_geometry), index_(index) { }
      /*!*******************************************************************************************
       * \brief Return data of the mapping in the tensor product of a 1-dimensional quadrature set.
       *
       * \tparam  npts  The number of evaluation points in a single direction
       * \tparam  T     The data type used for this operation
       ********************************************************************************************/
      template <std::size_t npts, typename T = double>
      Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, npts, T>
      mapping_tensor(const std::array<T, npts>& points_1d) const;
      /*!*******************************************************************************************
       * \brief   Return vertex of specified index of a hyperedge.
       ********************************************************************************************/
      Point<space_dimT> point(const unsigned int pt_index) const
      { 
        return hyGraph_geometry_.domain_info_.
                 points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][pt_index]];
      }
      /*!*******************************************************************************************
       * \brief   Return normal of specified index of a hyperedge.
       *
       * \todo    This function works only if hyperedge_dim == 1, so far.
       ********************************************************************************************/
      Point<space_dimT> normal(const unsigned int index) const
      {
        hy_assert( hyEdge_dimT == 1 , "This function has only been implemented for 1D edges." );
        hy_assert( index < 2 , "One dimensional edges can only have two normals." );

        Point<space_dimT> normal = 
          hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.
            domain_info_.points_hyEdge[index_][0]]
          - hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.
            domain_info_.points_hyEdge[index_][1]] ;
        normal /= norm_2(normal) * (2 * (index == 0) - 1);
        return normal;
      }
  }; // end of class hyEdge
  
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * \retval  hyEdge_dimT   The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of the space.
     *
     * \retval  space_dimT       The dimension of the space.
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
     * To receive a very general \c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The \c typedef implements the aspect, that a cubic hypergraph
     * topology is by default constructed by a std::vector that contains amounts of elements in the
     * different dimensions.
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
     * \brief   Get geometrical hyperedge of given index.
     *
     * This function returns the hyperedge of the given index, i.e., it returns the geometrical
     * hyperedge (\b not the topological information). The geometrical informatiom comprises the
     * indices of adjacent vertices (i.e. points) and information about their respective positions.
     *
     * \param   index       The index of the hyperedge to be returned.
     * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    const value_type operator[](const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
}; // end class File


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF Geometry::File
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// mapping_tensor
// -------------------------------------------------------------------------------------------------

template <unsigned int edim, unsigned int sdim, typename a , typename mapping_t, typename hyEdge_index_t>
template <std::size_t npts, typename T>
Tensor::MappingMultilinear<sdim, edim, npts, T>
File<edim, sdim, a, mapping_t,hyEdge_index_t>::hyEdge::mapping_tensor(const std::array<T, npts>& points_1d) const
{
  std::array<Point<sdim>, Hypercube<edim>::n_vertices()> vertices;
  for (unsigned int i=0;i<vertices.size();++i)  vertices[i] = point(i);
  return Tensor::MappingMultilinear<sdim, edim, npts, T>(vertices, points_1d);
}

} // end namespace Geometry
