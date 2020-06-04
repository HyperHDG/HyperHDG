#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/Hypercube.hxx>
#include <HyperHDG/DenseLA.hxx>
#include <HyperHDG/Topology/Cubic.hxx>

#include <HyperHDG/NodeDescriptor/Cubic.hxx>

#include <tensor_mapping.hxx>
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
 * \tparam  hyEdge_dimT    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for 
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dimT.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyEdge_index_t = unsigned int>
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
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
    private:
      /*!*******************************************************************************************
       * \brief   Points adjacent to the hyperedge.
       *
       * \todo    In the long run, do not store these.
       *
       * An array comprising the vertices (points) of a cubic hyperedge.
       ********************************************************************************************/
      std::array<Point<space_dimT>, Hypercube<hyEdge_dimT>::n_vertices()> points_;
    public:
      static constexpr unsigned int space_dim() { return space_dimT; }
      static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
      template <typename pt_coord_t>
      Point<space_dimT, pt_coord_t> map_ref_to_phys(const Point<hyEdge_dimT,pt_coord_t>& pt) const
      {
        Point<space_dimT> result = points_[0];
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          result += (float) pt[dim] * (points_[1<<dim] - points_[0]);
        return (Point<space_dimT,pt_coord_t>) result;
      }
      double area() const { return 1.; }
      Point<space_dimT> span_vec(const unsigned int index)
      {Point<space_dimT> a(1.); return a;}
      const SmallSquareMat<hyEdge_dimT> mat_r() { SmallSquareMat<hyEdge_dimT> a = diagonal<hyEdge_dimT,hyEdge_dimT,double>(1.); return a; }
      Point<hyEdge_dimT> local_normal(const unsigned int face)
      {
        Point<hyEdge_dimT> normal;
        normal[face/2] = 2. * (face % 2) - 1.;
        return normal;
      }
      double face_area(const unsigned int index) {return 1.;}
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
      hyEdge(const hyEdge_index_t index, const UnitCube& geometry)
      {
        
      };
    
      /*!*******************************************************************************************
       * \brief   Return vertex of specified index of a hyperedge.
       *
       * Return a \c Point describing the position of a vertex of a hyperedge.
       *
       * \retval  point           Point/Vertex of the hyperedge.
       ********************************************************************************************/
      Point<space_dimT> point(const unsigned int index) const  { return points_[index]; }

      /*!*******************************************************************************************
       * \brief Return data of the mapping in the tensor product of a one-dimensional quadrature set.
       *
       * \tparam npts: The number of evaluation points in a single direction
       * \tparam T: The data type used for this operation
       ********************************************************************************************/
      template <std::size_t npts, typename T = double>
      Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, npts, T>
      mapping_tensor(const std::array<T, npts>& points_1d) const;
    
      template<unsigned int n_sub_points, typename float_t>
      Point<space_dimT> lexicographic(unsigned int index, std::array<float_t, n_sub_points> points)
      {
        std::array<Point<space_dimT>, Hypercube<hyEdge_dimT>::n_vertices()> vertices;
        for (unsigned int i=0;i<vertices.size();++i)  vertices[i] = point(i);
    
        Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, n_sub_points, float_t>
          mapping(vertices, points);
        return mapping.lexicographic(index);
      }

      /*!*******************************************************************************************
       * \todo    Guido: If you have a clever idea for this, you can implement it. But this, I might
       *          also be able to do myself ;)
       ********************************************************************************************/
      Point<space_dimT> normal(const unsigned int index) const;

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
     * \brief   Number of elements per spatial dimension.
     *
     * A \c std::array comprising the number of elements in each spatial dimension.
     **********************************************************************************************/
    std::array<unsigned int, space_dimT> num_elements_;
    /*!*********************************************************************************************
     * \brief   Tensor product chain complex for elements.
     **********************************************************************************************/
    tpcc_chain<hyEdge_dimT, space_dimT, hyEdge_index_t> tpcc_;
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
     * \todo Check whether this is till a good idea. It hides important information behind a typedef.
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
     * \todo This is copied from the other constructor
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const constructor_value_type& num_elements)
    : tpcc_(create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements_))
    { 
      for (unsigned int dim = 0; dim < space_dimT; ++dim) num_elements_[dim] = num_elements[dim];
      tpcc_ = create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements_);
    }
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const Topology::Cubic<hyEdge_dimT,space_dimT>& other)
    : num_elements_(other.num_elements()),
      tpcc_(create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(other.num_elements()))
    { }
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
    { return hyEdge(index, *this); }
}; // end class UnitCube

  template <unsigned int edim, unsigned int sdim, typename hyEdge_index_t>
  template <std::size_t npts, typename T>
  Tensor::MappingMultilinear<sdim, edim, npts, T>
  UnitCube<edim, sdim, hyEdge_index_t>::hyEdge::mapping_tensor(const std::array<T, npts>& points_1d) const
  {
    std::array<Point<sdim>, Hypercube<edim>::n_vertices()> vertices;
    for (unsigned int i=0;i<vertices.size();++i)
      vertices[i] = point(i);

    return Tensor::MappingMultilinear<sdim, edim, npts, T>(vertices, points_1d);
  }

} // end namespace Geometry
