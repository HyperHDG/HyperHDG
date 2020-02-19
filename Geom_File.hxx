#ifndef GEOM_FILE_HXX
#define GEOM_FILE_HXX

#include <Point.hxx>
#include <tensor_mapping.hxx>
#include <array>


namespace Geometry
{

template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyEdge_index_t = unsigned int>
class File
{
  
  class hyEdge
  {
    public:
      static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
      static constexpr unsigned int space_dim() { return space_dimT; }
    private:
      const File& hyGraph_geometry_;
      const hyEdge_index_t index_;
    public:
      hyEdge ( const File& hyGraph_geometry, const hyEdge_index_t index )
      : hyGraph_geometry_(hyGraph_geometry), index_(index) { }

    /*!*******************************************************************************************
     * \brief Return data of the mapping in the tensor product of a one-dimensional quadrature set.
     *
     * \tparam npts: The number of evaluation points in a single direction
     * \tparam T: The data type used for this operation
     ********************************************************************************************/
      template <std::size_t npts, typename T = double>
      Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, npts, T>
      mapping_tensor(const std::array<T, npts>& points_1d) const;

      Point<space_dimT> point(const unsigned int pt_index) const
      { return hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][pt_index]]; }
      Point<space_dimT> normal(const unsigned int index) const
      {
        Point<space_dimT> normal;
        if (index == 0)  normal = hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][0]]
                                  - hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][1]] ;
        else             normal = hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][1]]
                                  - hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][0]] ;
        normal /= norm_2(normal);
        return normal;
      }
  }; // end of class hyEdge
  
  public:
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dimT; } // To be removed
    static constexpr unsigned int space_dimension() { return space_dimT; } // To be removed
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    static constexpr unsigned int space_dim() { return space_dimT; }
  private:
    const DomainInfo<hyEdge_dimT,space_dimT>& domain_info_;
  public:
    typedef hyEdge value_type;
    typedef Topology::File<hyEdge_dimT,space_dimT> constructor_value_type;
    File(const constructor_value_type& topology) : domain_info_(topology.domain_info()) { }
    const hyEdge get_hyEdge(const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
}; // end class File

  template <unsigned int edim, unsigned int sdim, typename hyEdge_index_t>
  template <std::size_t npts, typename T>
  Tensor::MappingMultilinear<sdim, edim, npts, T>
  File<edim, sdim, hyEdge_index_t>::hyEdge::mapping_tensor(const std::array<T, npts>& points_1d) const
  {
    std::array<Point<sdim>, Hypercube<edim>::n_vertices()> vertices;
    for (unsigned int i=0;i<vertices.size();++i)
      vertices[i] = point(i);

    return Tensor::MappingMultilinear<sdim, edim, npts, T>(vertices, points_1d);
  }

} // end namespace Geometry

#endif // end ifndef GEOM_FILE_HXX
