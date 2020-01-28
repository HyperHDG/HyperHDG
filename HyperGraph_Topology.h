#ifndef HYPERGRAPH_TOPOLOGY_H
#define HYPERGRAPH_TOPOLOGY_H

#include "TypeDefs.h"
#include "HyperEdge_Topology.h"
#include <array>
#include <vector>

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
 * @brief   Definition of the topology of a hypergraph --- Cubic HyperGraphs.
 *
 * \todo This is not what brief says. It is one special hypergraph.
 * \todo Do we use \\ or \@? Let us unify doxygen syntax
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
class HyperGraph_Cubic
{
  private:
    /*!*********************************************************************************************
     * @brief   Number of elements per spatial dimension.
     *
     * A @c std::array comprising the number of elements in each spatial dimension.
     **********************************************************************************************/
    std::array<unsigned int, space_dim> num_elements_;
    /*!*********************************************************************************************
     * @brief   Total amount of hyperedges.
     *
     * The number of hyperedges that form the hypergraph. This information is needed to allow to go
     * through all hyperedges and execute some code. The number of hyperedges can be computed from
     * the @c std::array @c num_elements_.
     **********************************************************************************************/
    hyperedge_index_type n_hyperedges_;
    /*!*********************************************************************************************
     * @brief   Total amount of hypernodes.
     *
     * The number of hypernodes that make up the hypergraph. This information is needed to have the
     * appropriate version of a @c HyperNodeFactory. It can be vomputed from the \c std::array
     * @c num_elements_.
     **********************************************************************************************/
    hypernode_index_type n_hypernodes_;
  public:
    /*!*********************************************************************************************
     * @brief   Defines the return value of the class.
     *
     * The @c class @c HyperGraph_Cubic defines the topology of the hypergraph. It "contains" the
     * different hyperedges (that actually are constructed everytime access is needed from e.g. the
     * solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef HyperEdge_Cubic<hyperedge_dim, space_dim> value_type;
    /*!*********************************************************************************************
     * @brief   Defines the value type of input argument for standard constructor.
     *
     * To receive a very general @c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The @c typedef implements the aspect, that a cubic hypergraph
     * topology is by default constructed by a std::vector that contains amounts of elements in the
     * different dimensions.
     **********************************************************************************************/
    typedef std::vector<int> constructor_value_type;
    /*!*********************************************************************************************
     * @brief   Construct a cubic hypergraph from a @c std::vector.
     *
     * Constructs a hypergraph from a @c std::vector containing the elementens per spatial dimension
     * which is given as input data. If the input vector is shorter that @c space_dim, the remaining
     * amounts of elemnts are assumed to be equal to zero. If the vector is longer than
     * @c space_dim, the first @c space_dim entries are considered only.
     * 
     * @todo    If the vector is too short, an error is thrown in the test program and the behavior
     *          is undefined for Python (most likely an error is thrown, too) at the moment.
     *
     * @param   num_elements    A @c std::vector containing number of elements per dimension.
     **********************************************************************************************/
    HyperGraph_Cubic(const constructor_value_type& num_elements);
    /*!*********************************************************************************************
     * @brief   Construct a cubic hypergraph from a @c std::array.
     *
     * Constructs a hypergraph from a @c std::array containing the elementens per spatial dimension
     * which is given as input data. The array has the correct length (as ensured by the involved
     * template parametzer @c space_dim.
     *
     * @param   num_elements    A @c std::array containing number of elements per spatial dimension.
     **********************************************************************************************/
    HyperGraph_Cubic(const std::array<unsigned int, space_dim>& num_elements);
    /*!*********************************************************************************************
     * @brief   Construct a hypergraph from another hypergraph.
     *
     * Create a (value based) copy of another hypergraph.
     *
     * @param   other           Hypergraph to be copied.
     **********************************************************************************************/
    HyperGraph_Cubic(const HyperGraph_Cubic<hyperedge_dim,space_dim>& other);
    /*!*********************************************************************************************
     * @brief   Get topological hyperedge of given index.
     *
     * This function returns the hyperedge of the given index, i.e., it returns the topological
     * hyperedge (@b not the geometrical information). The topological informatiom comprises the
     * indices of adjacent hypernodes and information about their respective orientations.
     *
     * @param   index           The index of the hyperedge to be returned.
     * @retval  hyperedge       Topological information on the hyperedge (cf. @c value_type).
     **********************************************************************************************/
    const HyperEdge_Cubic<hyperedge_dim, space_dim>
      get_hyperedge(const hyperedge_index_type index) const;
    /*!*********************************************************************************************
     * @brief   Read the array of elements per dimensions.
     *
     * @retval  num_elements    A @c std::array containing the elements in the repective dimension.
     **********************************************************************************************/
    const std::array<unsigned int, space_dim>& num_elements() const;
    /*!*********************************************************************************************
     * @brief   Returns the number of hyperedges making up the hypergraph.
     *
     * @retval  n_hyperedges    The total amount of hyperedges of a hypergraph.
     **********************************************************************************************/
    const hyperedge_index_type n_hyperedges() const;
    /*!*********************************************************************************************
     * @brief   Returns the number of hypernodes making up the hypergraph.
     *
     * @retval  n_hypernodes    The total amount of hypernodes of a hypergraph.
     **********************************************************************************************/
    const hypernode_index_type n_hypernodes() const;
    
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
}; // end of class HyperGraph_Cubic

} // end of namespace Topology

#endif // end of ifndef HYPERGRAPH_TOPOLOGY_H