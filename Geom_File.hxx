#ifndef GEOM_FILE_HXX
#define GEOM_FILE_HXX

#include <TypeDefs.hxx>
#include <Point.hxx>
#include <Topo_File.hxx>
#include <array>


namespace Geometry
{

template <unsigned int hyEdge_dim, unsigned int space_dim>
class File
{
  
  class hyEdge
  {
    private:
      const File& hyGraph_geometry_;
      const hyEdge_index_t index_;
    public:
      hyEdge ( const File& hyGraph_geometry, const hyEdge_index_t index )
      : hyGraph_geometry_(hyGraph_geometry), index_(index) { }
                               
      Point<space_dim> point(const unsigned int pt_index) const
      { return hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][pt_index]]; }
      Point<space_dim> normal(const unsigned int index) const
      {
        Point<space_dim> normal;
        if (index == 0)  normal = hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][0]]
                                  - hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][1]] ;
        else             normal = hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][1]]
                                  - hyGraph_geometry_.domain_info_.points[hyGraph_geometry_.domain_info_.points_hyEdge[index_][0]] ;
        normal /= norm_2(normal);
        return normal;
      }
  }; // end of class hyEdge
  
  private:
    const DomainInfo<hyEdge_dim,space_dim>& domain_info_;
  public:
    typedef hyEdge value_type;
    typedef Topology::File<hyEdge_dim,space_dim> constructor_value_type;
    File(const constructor_value_type& topology) : domain_info_(topology.domain_info()) { }
    const hyEdge get_hyEdge(const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
    
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr unsigned int space_dimension() { return space_dim; }
}; // end class File

} // end namespace Geometry

#endif // end ifndef GEOM_FILE_HXX
