#ifndef HYPEREDGE_TOPOLOGY_H
#define HYPEREDGE_TOPOLOGY_H

#include "TypeDefs.h"
#include <array>

/*!*************************************************************************************************
 * @brief   A namespace containing different classes describing hypergraph topologies.
 *
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain @f$\Omega@f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/
namespace Topology
{

/*!*************************************************************************************************
 * @brief   Definition of the topology of a hypergraph's edges --- Cubic HyperGraph's edges.
 * 
 * @todo    Both private arrays are filled when the hyperedge is constructed. Lazy evaluation might
 *          be an important aspect here. What do you think?
 * 
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain @f$\Omega@f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * @tparam  hyperedge_dim   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * @tparam  space_dim       The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyperedge_dim.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/
template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperEdge_Cubic
{
  private:
    /*!*********************************************************************************************
     * @brief   Indices of the hypernodes adjacent to the hyperedge.
     *
     * A @c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
     **********************************************************************************************/
    std::array<hypernode_index_type, 2*hyperedge_dim> hypernode_indices_;
    /*!*********************************************************************************************
     * @brief   Orientation of the hypernode.
     * 
     * @todo    Do we want to change this (cf. detailed description)? This array also does not have
     *          a getter function!
     * 
     * A @c std::array comprising the orientation of each hypernode. In HyperGraph_Cubic, all edges
     * are assumed to have the correct orientation and this array is irrelevant. However, this is
     * possible to change for different applications.
     **********************************************************************************************/
    std::array<unsigned int, 2*hyperedge_dim> correct_hypernode_orientation_;
  public:
    /*!*********************************************************************************************
     * @brief   Construct a cubic hyperedge from its index and a \c std::array of elements in each
     *          spatial dimension.
     *
     * Constructs a hyperedge from a @c std::array containing the elementens per spatial dimension
     * which is given as input data and the index of the hyperedge to be constructed.
     * 
     * @param   index           The index of the hyperedge to be created.
     * @param   num_elements    A @c std::array containing number of elements per dimension.
     **********************************************************************************************/
    HyperEdge_Cubic(const hyperedge_index_type index,
                    const std::array<unsigned int, space_dim>& num_elements);
    /*!*********************************************************************************************
     * @brief   Return indices of hypernodes adjacent to the hyperedge.
     *
     * Return a @c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
     *
     * @retval  hypernode_indeices  Topological information on the hyperedge (cf. @c value_type).
     **********************************************************************************************/
    const std::array<hypernode_index_type, 2*hyperedge_dim>& get_hypernode_indices() const;
}; // end of class HyperEdge_Cubic

} // end of namespace Topology

#endif // end of ifndef HYPEREDGE_TOPOLOGY_H
