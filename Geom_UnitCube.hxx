/*!*************************************************************************************************
 * \file    Geom_UnitCube.hxx
 * \brief   Define geometry of hypergraphs that form unit cubes.
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
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/

#ifndef GEOM_UNITCUBE_HXX
#define GEOM_UNITCUBE_HXX

#include <Hypercube.hxx>
#include <TypeDefs.hxx>
#include <Point.hxx>
#include <Topo_Cubic.hxx>
#include <array>

/*!*************************************************************************************************
 * \brief   A namespace containing different classes describing hypergraph geometries.
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
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
namespace Geometry
{

/*!*************************************************************************************************
 * \brief   Definition of the topology of a hypergraph --- Cubic HyperGraphs that are unit cubes.
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
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for 
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dim     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dim.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim>
class UnitCube
{
  
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's hyperedge --- Edges of cubic HyperGraphs
   *          that form unit cubes.
   * 
   * \todo    We do not use this information for our simulations, but for plotting only. Thus, the
   *          results do not really show the truth. A HyperGraph consisiting of HyperEdges that are
   *          unit cubes is calculated and a HyperGraph that is a cube and made up of smaller
   *          quadrilaterals (in 2D) is plotted!
   * 
   * \todo    Functions returning the Jacobian / its det, etc are not yet implemented. Should they
   *          be implemented in the future or should their interfaces (commented at the moment) be
   *          removed from this file?
   * 
   * \todo    The points/vertices are computed when the HyperEdge is constructed. Lazy evaluation
   *          might be a relevant aspect here. What do you think?
   * 
   * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
   **************************************************************************************************/
  class hyEdge
  {
    private:
      /*!*******************************************************************************************
       * \brief   Points adjacent to the hyperedge.
       *
       * A \c std::array comprising the vertices (points) of a cubic hyperedge.
       ********************************************************************************************/
      std::array<Point<space_dim>, Hypercube<hyEdge_dim>::n_vertices()> points_;
    public:
      /*!*******************************************************************************************
       * \brief   Construct a cubic hyperedge from its index and a \c std::array of elements in each
       *          spatial dimension.
       *
       * \todo    Guido: Please, implement function that constructs the hyperedge of a given index.
       *          A prototype of this function is located in the cxx file, where you could also
       *          insert the new function.
       * 
       * Constructs a hyperedge from a \c std::array containing the elementens per spatial dimension
       * which is given as input data and the index of the hyperedge to be constructed.
       * 
       * \param   index           The index of the hyperedge to be created.
       * \param   num_elements    A \c std::array containing number of elements per dimension.
       ********************************************************************************************/
      hyEdge(const hyEdge_index_t index, const std::array<unsigned int, space_dim>& num_elements);
      /*!*******************************************************************************************
       * \brief   Return vertex of specified index of a hyperedge.
       *
       * Return a \c Point describing the position of a vertex of a hyperedge.
       *
       * \retval  point           Point/Vertex of the hyperedge.
       ********************************************************************************************/
      Point<space_dim> point(const unsigned int index) const
      { return points_[index]; }
      
      /*!*******************************************************************************************
       * \todo    Guido: If you have a clever idea for this, you can implement it. But this, I might
       *          also be able to do myself ;)
       ********************************************************************************************/
      Point<space_dim> normal(const unsigned int index) const;

//    std::vector<double> abs_det_of_jacobian_at_quad
//      (const std::vector<double>& local_quadrature) const;
//    std::vector< std::vector<double> > inv_of_transposed_jacobian_at_quad
//      (const std::vector<double>& local_quadrature) const;
  }; // end of class hyEdge
  
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * \retval  hyEdge_dim   The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of the space.
     *
     * \retval  space_dim       The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dimension() { return space_dim; }
  private:
    /*!*********************************************************************************************
     * \brief   Number of elements per spatial dimension.
     *
     * A \c std::array comprising the number of elements in each spatial dimension.
     **********************************************************************************************/
    std::array<unsigned int, space_dim> num_elements_;
  public:
    /*!*********************************************************************************************
     * \brief   Defines the return value of the class.
     *
     * The \c class \c UnitCube defines the geometry of the hypergraph. It contains
     * the different hyperedges (that actually are constructed everytime access is needed from e.g.
     * the solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef hyEdge value_type;
    /*!*********************************************************************************************
     * \brief   Defines the value type of input argument for standard constructor.
     *
     * To receive a very general \c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The \c typedef implements the aspect, that a cubic hypergraph
     * topology whic is a unit cube is by default constructed by a std::vector that contains amounts
     * of elements in the different dimensions.
     **********************************************************************************************/
    typedef std::vector<unsigned int> constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const constructor_value_type& num_elements)
    { for (unsigned int dim = 0; dim < space_dim; ++dim) num_elements_[dim] = num_elements[dim]; }
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const Topology::Cubic<hyEdge_dim,space_dim>& other)
    : num_elements_(other.num_elements()) { }
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
    const hyEdge get_hyEdge(const hyEdge_index_t index) const
    { return hyEdge(index, num_elements_); }
}; // end class UnitCube

} // end namespace Geometry

#endif // end ifndef GEOM_UNITCUBE_HXX
