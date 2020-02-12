#ifndef TOPO_FILE_HXX
#define TOPO_FILE_HXX

#include "TypeDefs.hxx"
#include "ReadDomain.hxx"
#include <array>
#include <string>


namespace Topology
{  

template <unsigned int hyEdge_dimT, unsigned int space_dimT>
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
      template < typename hyNode_index_t = unsigned int >
      const std::array<hyNode_index_t, 2*hyEdge_dimT>& get_hyNode_indices() const
      { return hyGraph_topology_.domain_info_.hyNodes_hyEdge[index_]; }
  }; // end of class hyEdge
  
  public:
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    static constexpr unsigned int space_dim() { return space_dimT; }
  
  private:
    const DomainInfo<hyEdge_dimT,space_dimT> domain_info_;
    
  public:
    typedef hyEdge value_type;
    typedef std::string& constructor_value_type;
    
    File(const constructor_value_type& filename)
    : domain_info_(read_domain<hyEdge_dimT,space_dimT>(filename)) { }
    File(const File<hyEdge_dimT,space_dimT>& other)
    : domain_info_(other.domain_info) { }
    const value_type operator[](const hyEdge_index_t index) const
    {
      hy_assert( index < domain_info_.n_hyEdges && index >= 0 ,
                 "Index must be non-negative and smaller than " << domain_info_.n_hyEdges <<
                 " (which is the amount of hyperedges). It was " << index << "!" );
      return hyEdge(*this, index);
    }
    const hyEdge_index_t n_hyEdges() const { return domain_info_.n_hyEdges; }
    template < typename hyNode_index_t = unsigned int >
    const hyNode_index_t n_hyNodes() const { return domain_info_.n_hyNodes; }
    const DomainInfo<hyEdge_dimT,space_dimT>& domain_info() const { return domain_info_; }
}; // end of class File

} // end of namespace Topology

#endif // end of ifndef TOPO_FILE_HXX
