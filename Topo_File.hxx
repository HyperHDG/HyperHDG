#ifndef TOPO_FILE_HXX
#define TOPO_FILE_HXX

#include "TypeDefs.hxx"
#include "ReadDomain.hxx"
#include <array>
#include <string>


namespace Topology
{  

template <unsigned int hyEdge_dim, unsigned int space_dim>
class File
{
  class hyEdge
  {
    private:
      const File& hyGraph_topology_;
      const hyEdge_index_t index_;
    public:
      hyEdge ( const File& hyGraph_topology, const hyEdge_index_t index )
      : hyGraph_topology_(hyGraph_topology), index_(index) { }
      const std::array<hyNode_index_t, 2*hyEdge_dim>& get_hyNode_indices() const
      { return hyGraph_topology_.domain_info_.hyNodes_hyEdge[index_]; }
  }; // end of class hyEdge
  
  private:
    const DomainInfo<hyEdge_dim,space_dim> domain_info_;
    
  public:
    typedef hyEdge value_type;
    typedef std::string& constructor_value_type;
    
    File(const constructor_value_type& filename)
    : domain_info_(read_domain<hyEdge_dim,space_dim>(filename)) { }
    File(const File<hyEdge_dim,space_dim>& other)
    : domain_info_(other.domain_info) { }
    const hyEdge get_hyEdge(const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
    const hyEdge_index_t n_hyEdges() const { return domain_info_.n_hyEdges; }
    const hyNode_index_t n_hyNodes() const { return domain_info_.n_hyNodes; }
    const DomainInfo<hyEdge_dim,space_dim>& domain_info() const { return domain_info_; }
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; };
    static constexpr unsigned int space_dimension() { return space_dim; };
}; // end of class File

} // end of namespace Topology

#endif // end of ifndef TOPO_FILE_HXX
