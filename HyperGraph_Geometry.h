#ifndef HYPERGRAPH_GEOMETRY_H
#define HYPERGRAPH_GEOMETRY_H

#include "TypeDefs.h"
#include "HyperEdge_Geometry.h"
#include "Topo_Cubic.h"
#include <array>

/*!*************************************************************************************************
 * @brief   A namespace containing different classes describing hypergraph geometries.
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
namespace Geometry
{

/*!*************************************************************************************************
 * @brief   Definition of the topology of a hypergraph --- Cubic HyperGraphs that are unit cubes.
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
class HyperGraph_Cubic_UnitCube
{
  private:
    /*!*********************************************************************************************
     * @brief   Number of elements per spatial dimension.
     *
     * A \c std::array comprising the number of elements in each spatial dimension.
     **********************************************************************************************/
    std::array<unsigned int, space_dim> num_elements_;
  public:
    /*!*********************************************************************************************
     * @brief   Defines the return value of the class.
     *
     * The @c class @c HyperGraph_Cubic_UnitCube defines the geometry of the hypergraph. It contains
     * the different hyperedges (that actually are constructed everytime access is needed from e.g.
     * the solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> value_type;
    /*!*********************************************************************************************
     * @brief   Defines the value type of input argument for standard constructor.
     *
     * To receive a very general @c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The @c typedef implements the aspect, that a cubic hypergraph
     * topology whic is a unit cube is by default constructed by a std::vector that contains amounts
     * of elements in the different dimensions.
     **********************************************************************************************/
    typedef std::vector<int> constructor_value_type;
    /*!*********************************************************************************************
     * @brief   Construct a cubic that describes a cube hypergraph from a @c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a @c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * @param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    HyperGraph_Cubic_UnitCube(const constructor_value_type& num_elements);
    /*!*********************************************************************************************
     * @brief   Construct a cubic that describes a cube hypergraph from a @c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a @c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * @param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    HyperGraph_Cubic_UnitCube(const Topology::Cubic<hyperedge_dim,space_dim>& other);
    /*!*********************************************************************************************
     * @brief   Get geometrical hyperedge of given index.
     *
     * This function returns the hyperedge of the given index, i.e., it returns the geometrical
     * hyperedge (@b not the topological information). The geometrical informatiom comprises the
     * indices of adjacent vertices (i.e. points) and information about their respective positions.
     *
     * @param   index       The index of the hyperedge to be returned.
     * @retval  hyperedge   Geometrical information on the hyperedge (cf. @c value_type).
     **********************************************************************************************/
    const HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> get_hyperedge
      (const hyperedge_index_type index) const;
      
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * @retval  hyperedge_dim   The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyperedge_dimension() { return hyperedge_dim; };
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the dimension of the space.
     *
     * @retval  space_dim       The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dimension() { return space_dim; };
}; // end class HyperGraph_Cubic_UnitCube

} // end namespace Geometry

#endif // end ifndef HYPERGRAPH_GEOMETRY_H
