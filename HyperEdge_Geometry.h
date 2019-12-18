#ifndef HYPEREDGE_GEOMETRY_H
#define HYPEREDGE_GEOMETRY_H

#include "TypeDefs.h"
#include "Point.h"
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
 * @brief   Definition of the topology of a hypergraph's hyperedge --- Edges of cubic HyperGraphs
 *          that form unit cubes.
 * 
 * @todo    We do not use this information for our simulations, but for plotting only. Thus, the
 *          results do not really show the truth. A HyperGraph consisiting of HyperEdges that are
 *          unit cubes is calculated and a HyperGraph that is a cube and made up of smaller
 *          quadrilaterals (in 2D) is plotted!
 * 
 * @todo    Functions returning the Jacobian / its det, etc are not yet implemented. Should they be
 *          implemented in the future or should their interfaces (commented at the moment) be
 *          removed from this file?
 * 
 * @todo    The points/vertices are computed when the HyperEdge is constructed. Lazy evaluation
 *          might be a relevant aspect here. What do you think?
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
class HyperEdge_Cubic_UnitCube
{
  private:
    /*!*********************************************************************************************
     * @brief   Points adjacent to the hyperedge.
     *
     * A \c std::array comprising the vertices (points) of a cubic hyperedge.
     **********************************************************************************************/
    std::array<Point<space_dim>, 2*hyperedge_dim> points_;
  public:
    /*!*********************************************************************************************
     * @brief   Construct a cubic hyperedge from its index and a \c std::array of elements in each
     *          spatial dimension.
     *
     * Constructs a hyperedge from a \c std::array containing the elementens per spatial dimension
     * which is given as input data and the index of the hyperedge to be constructed.
     * 
     * @param   index           The index of the hyperedge to be created.
     * @param   num_elements    A @c std::array containing number of elements per dimension.
     **********************************************************************************************/
    HyperEdge_Cubic_UnitCube(const hyperedge_index_type index,
                             const std::array<unsigned int, space_dim>& num_elements);
    /*!*********************************************************************************************
     * @brief   Return vertex of specified index of a hyperedge.
     *
     * Return a \c Point describing the position of a vertex of a hyperedge.
     *
     * @retval  point           Point/Vertex of the hyperedge.
     **********************************************************************************************/
    Point<space_dim> point(unsigned int index) const;

//    std::vector<double> abs_det_of_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
//    std::vector< std::vector<double> > inv_of_transposed_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
}; // end of class HyperEdge_Cubic_UnitCube

} // end of namespace Geometry

#endif // end of ifndef HYPEREDGE_GEOMETRY_H
