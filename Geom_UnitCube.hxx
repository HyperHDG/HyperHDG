#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/Hypercube.hxx>
#include <HyperHDG/DenseLA.hxx>
#include <Topo_Cubic.hxx>
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
       * \todo In the long run, do not store these.
       *
       * An array comprising the vertices (points) of a cubic hyperedge.
       ********************************************************************************************/
      std::array<Point<space_dimT>, Hypercube<hyEdge_dimT>::n_vertices()> points_;
    public:
      static constexpr unsigned int space_dim() { return space_dimT; }
      static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
      Point<space_dimT> map_ref_to_phys(const Point<hyEdge_dimT>& pt) const
      {Point<space_dimT> a; return a;}
      double area() const { return 1.; }
      Point<space_dimT> span_vec(const unsigned int index)
      {Point<space_dimT> a(1.); return a;}
      const SmallSquareMat<hyEdge_dimT> mat_r() { SmallSquareMat<hyEdge_dimT> a = diagonal<hyEdge_dimT,hyEdge_dimT,double>(1.); return a; }
      Point<hyEdge_dimT> hyEdge_dim_normal(const unsigned int face)
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
      hyEdge(const hyEdge_index_t index, const std::array<unsigned int, space_dimT>& num_elements);
    
      /*!*******************************************************************************************
       * \brief   Return vertex of specified index of a hyperedge.
       *
       * Return a \c Point describing the position of a vertex of a hyperedge.
       *
       * \retval  point           Point/Vertex of the hyperedge.
       ********************************************************************************************/
      Point<space_dimT> point(const unsigned int index) const
      { return points_[index]; }

      /*!*******************************************************************************************
       * \brief Return data of the mapping in the tensor product of a one-dimensional quadrature set.
       *
       * \tparam npts: The number of evaluation points in a single direction
       * \tparam T: The data type used for this operation
       ********************************************************************************************/
    template <std::size_t npts, typename T = double>
      Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, npts, T>
      mapping_tensor(const std::array<T, npts>& points_1d) const;
    
      template<unsigned int n_sub_points>
      Point<space_dimT> lexicographic(unsigned int index)
      {
        std::array<Point<space_dimT>, Hypercube<hyEdge_dimT>::n_vertices()> vertices;
        for (unsigned int i=0;i<vertices.size();++i)  vertices[i] = point(i);
    
        std::array<double,n_sub_points> helper;
        Tensor::MappingMultilinear<space_dimT, hyEdge_dimT, n_sub_points, double>
          mapping(vertices, helper);
        return mapping.lexicographic(index);
      }

      /*!*******************************************************************************************
       * \todo    Guido: If you have a clever idea for this, you can implement it. But this, I might
       *          also be able to do myself ;)
       ********************************************************************************************/
      Point<space_dimT> normal(const unsigned int index) const;

//    std::vector<double> abs_det_of_jacobian_at_quad
//      (const std::vector<double>& local_quadrature) const;
//    std::vector< std::vector<double> > inv_of_transposed_jacobian_at_quad
//      (const std::vector<double>& local_quadrature) const;
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
    { for (unsigned int dim = 0; dim < space_dimT; ++dim) num_elements_[dim] = num_elements[dim]; }
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const Topology::Cubic<hyEdge_dimT,space_dimT>& other)
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
    const value_type operator[](const hyEdge_index_t index) const
    { return hyEdge(index, num_elements_); }
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






/*
 * HyperEdge functions!
 */


template
< 
  unsigned int space_dimT, typename hyEdge_index_t = unsigned int, 
  typename hyNode_index_t = hyEdge_index_t
>
std::array<Point<space_dimT>, 2> line_to_points(const std::array<unsigned int, space_dimT>& num_lines, const hyEdge_index_t index)
{
  hy_assert( num_lines.size() == space_dimT , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  std::array<Point<space_dimT>, 2> point_indices;
  point_indices.fill(Point<space_dimT>());
  
  std::array<hyEdge_index_t, space_dimT> num_lines_with_orientation;
  num_lines_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
      if (dim_m == dim_n)  num_lines_with_orientation[dim_m] *= num_lines[dim_n];
      else                 num_lines_with_orientation[dim_m] *= num_lines[dim_n] + 1;
  
  for ( orientation = 0;
        number_with_lower_orientation + num_lines_with_orientation[orientation] <= index ;
        ++orientation)
  {
    number_with_lower_orientation += num_lines_with_orientation[orientation];
    hy_assert( orientation <= space_dimT , "Orientation is a space_dimT and connot exceed it." );
  }
  
  std::array<hyEdge_index_t, space_dimT> local_indices;
  local_indices.fill(0);
  
  index_helper -= number_with_lower_orientation;
  for (unsigned int dim = 0; dim < space_dimT; ++dim)
  {
    num_elements_in_direction = num_lines[(dim + orientation) % space_dimT] + (dim != 0);
    local_indices[(dim + orientation) % space_dimT] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dimT];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
  for(unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
  {
    double helper = num_lines[dim_m];
    point_indices[0][dim_m] = local_indices[dim_m] / helper;
    if (dim_m == orientation)  point_indices[1][dim_m] = (local_indices[dim_m] + 1) / helper;
    else                       point_indices[1][dim_m] = local_indices[dim_m] / helper;
  }
  
  return point_indices;
}


template
< 
  unsigned int space_dimT, typename hyEdge_index_t = unsigned int, 
  typename hyNode_index_t = hyEdge_index_t
>
std::array<hyNode_index_t, 4> square_to_line_index(const std::array<unsigned int, space_dimT>& num_squares, const hyEdge_index_t index)
{
  hy_assert( num_squares.size() == space_dimT , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  std::array<hyNode_index_t, 4> line_indices;
  line_indices.fill(0);
  
  std::array<hyEdge_index_t, space_dimT> num_squares_with_orientation;
  num_squares_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
      if ( dim_m != dim_n )     num_squares_with_orientation[dim_m] *= num_squares[dim_n];
      else if (space_dimT == 2)  num_squares_with_orientation[dim_m] *= num_squares[dim_n];
      else if (space_dimT == 3)  num_squares_with_orientation[dim_m] *= num_squares[dim_n] + 1;
  
  for ( orientation = 0;
        number_with_lower_orientation + num_squares_with_orientation[orientation] <= index ;
        ++orientation)
  {
    number_with_lower_orientation += num_squares_with_orientation[orientation];
    hy_assert( orientation <= space_dimT , "Orientation is a space_dimT and connot exceed it." );
  }
  
  std::array<hyEdge_index_t, space_dimT> local_indices;
  local_indices.fill(0);
  
  index_helper -= number_with_lower_orientation;
  for (unsigned int dim = 0; dim < space_dimT; ++dim)
  {
    if (space_dimT == 3)  num_elements_in_direction = num_squares[(dim + orientation) % space_dimT] + (dim == 0);
    else                 num_elements_in_direction = num_squares[(dim + orientation) % space_dimT];
    local_indices[(dim + orientation) % space_dimT] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dimT];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
  std::array<hyEdge_index_t, space_dimT> num_lines_with_orientation;
  num_lines_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
      if (dim_m == dim_n)  num_lines_with_orientation[dim_m] *= num_squares[dim_n];
      else                 num_lines_with_orientation[dim_m] *= num_squares[dim_n] + 1;
  
  unsigned int local_line_index = 0;
  for (unsigned int line_orientation = 0; line_orientation < space_dimT; ++line_orientation)
  {
    if (space_dimT == 3 && line_orientation == orientation)  continue;
    number_with_lower_orientation = 0;
    for (unsigned int dim = 0; dim < line_orientation; ++dim)
      number_with_lower_orientation += num_lines_with_orientation[dim];
    for (int dim = space_dimT - 1; dim >= 0; --dim)
    {
      num_elements_in_direction = num_squares[(dim + line_orientation) % space_dimT] + (dim != 0);
      line_indices[local_line_index] *= num_elements_in_direction;
      line_indices[local_line_index] += local_indices[(dim + line_orientation) % space_dimT];
      line_indices[local_line_index + 1] *= num_elements_in_direction;
      if (space_dimT == 2) line_indices[local_line_index + 1] += local_indices[(dim + line_orientation) % space_dimT] + (dim != 0);
      else  line_indices[local_line_index + 1] += local_indices[(dim + line_orientation) % space_dimT] 
                                                  + ((dim + line_orientation) % space_dimT != orientation && dim != 0);
    }
    line_indices[local_line_index] += number_with_lower_orientation;
    line_indices[local_line_index + 1] += number_with_lower_orientation;
    local_line_index += 2;
  }
  
  return line_indices;
}


template
< 
  unsigned int space_dimT, typename hyEdge_index_t = unsigned int, 
  typename hyNode_index_t = hyEdge_index_t
>
std::array<hyNode_index_t, 6> cube_to_square_index(const std::array<unsigned int, space_dimT>& num_cubes, const hyEdge_index_t index)
{
  hy_assert( num_cubes.size() == space_dimT , "The size of the handed over parmeter does not fit!" );
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  std::array<hyNode_index_t, 6> square_indices;
  square_indices.fill(0);
  
  std::array<hyEdge_index_t, space_dimT> local_indices;
  local_indices.fill(0);
  
  for (unsigned int dim = 0; dim < space_dimT; ++dim)
  {
    num_elements_in_direction = num_cubes[dim];
    local_indices[dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[dim];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
  std::array<hyEdge_index_t, space_dimT> num_squares_with_orientation;
  num_squares_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
      if ( dim_m != dim_n )     num_squares_with_orientation[dim_m] *= num_cubes[dim_n];
      else if (space_dimT == 2)  num_squares_with_orientation[dim_m] *= num_cubes[dim_n];
      else if (space_dimT == 3)  num_squares_with_orientation[dim_m] *= num_cubes[dim_n] + 1;
  
  unsigned int local_square_index = 0;
  for (unsigned int square_orientation = 0; square_orientation < space_dimT; ++square_orientation)
  {
    number_with_lower_orientation = 0;
    for (unsigned int dim = 0; dim < square_orientation; ++dim)
      number_with_lower_orientation += num_squares_with_orientation[dim];
    for (int dim = space_dimT - 1; dim >= 0; --dim)
    {
      num_elements_in_direction = num_cubes[(dim + square_orientation) % space_dimT] + (dim == 0);
      square_indices[local_square_index] *= num_elements_in_direction;
      square_indices[local_square_index] += local_indices[(dim + square_orientation) % space_dimT];
      square_indices[local_square_index + 1] *= num_elements_in_direction;
      square_indices[local_square_index + 1] += local_indices[(dim + square_orientation) % space_dimT] + (dim == 0);
    }
    square_indices[local_square_index] += number_with_lower_orientation;
    square_indices[local_square_index + 1] += number_with_lower_orientation;
    local_square_index += 2;
  }
  
  return square_indices;
}

template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyE>
UnitCube<hyEdge_dimT,space_dimT, hyE>::hyEdge::
hyEdge(const hyE index, const std::array<unsigned int, space_dimT>& num_elements)
{
  if constexpr ( hyEdge_dimT == 1 )
  {
    points_ = line_to_points<space_dimT>(num_elements, index);
  }
  else if constexpr ( hyEdge_dimT == 2 )
  {
//    hy_assert( 0 == 1 , "This should not be executed since elasticity is not yet defined for real HYPERgraphs!" );
    std::array<hyE, 4> line_indices = square_to_line_index<space_dimT>(num_elements, index);
    std::array<Point<space_dimT>, 2> points1 = line_to_points<space_dimT>(num_elements, line_indices[0]);
    std::array<Point<space_dimT>, 2> points2 = line_to_points<space_dimT>(num_elements, line_indices[1]);
    points_[0] = points1[0];  points_[1] = points1[1];  points_[2] = points2[0];  points_[3] = points2[1];
  }
  else if constexpr ( hyEdge_dimT == 3 )
  {
//    hy_assert( 0 == 1 , "This should not be executed since elasticity is not yet defined for real HYPERgraphs!" );
    std::array<hyE, 6> square_indices = cube_to_square_index<space_dimT>(num_elements, index);
    std::array<hyE, 4> line_indices1 = square_to_line_index<space_dimT>(num_elements, square_indices[0]);
    std::array<hyE, 4> line_indices2 = square_to_line_index<space_dimT>(num_elements, square_indices[1]);
    std::array<Point<space_dimT>, 2> points1 = line_to_points<space_dimT>(num_elements, line_indices1[0]);
    std::array<Point<space_dimT>, 2> points2 = line_to_points<space_dimT>(num_elements, line_indices1[1]);
    std::array<Point<space_dimT>, 2> points3 = line_to_points<space_dimT>(num_elements, line_indices2[0]);
    std::array<Point<space_dimT>, 2> points4 = line_to_points<space_dimT>(num_elements, line_indices2[1]);
    points_[0] = points1[0];  points_[1] = points1[1];  points_[2] = points2[0];  points_[3] = points2[1];
    points_[4] = points3[0];  points_[5] = points3[1];  points_[6] = points4[0];  points_[7] = points4[1];
  }
//  sort(points_.begin(), points_.end());
}


template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyE>
Point<space_dimT> UnitCube<hyEdge_dimT,space_dimT,hyE>::hyEdge::
normal(unsigned int index) const
{
  Point<space_dimT> normal;
  if (index == 0)  normal = points_[0] - points_[1];
  else             normal = points_[1] - points_[0];
  normal /= norm_2(normal);
  return normal;
}

/*
const std::array<hyNode_index_t, 2*hyEdge_dimT>&
UnitCube<hyEdge_dimT,space_dimT>::get_hyNode_indices() const
{
  return hyNode_indices_;
}


template <unsigned int hyEdge_dimT, unsigned int space_dimT>
std::vector<double> UnitCube<hyEdge_dimT,space_dimT>::
abs_det_of_jacobian_at_quad(const vector<double>& local_quadrature) const
{
  vector<double> det_at_quad(local_quadrature.size(), 1.);
  return det_at_quad;
}


template <unsigned int hyEdge_dimT, unsigned int space_dimT>
vector< vector<double> > UnitCube<hyEdge_dimT,space_dimT>::
inv_of_transposed_jacobian_at_quad(const vector<double>& local_quadrature) const
{
  vector< vector<double> > jac_at_quad(local_quadrature.size());
  for_each(jac_at_quad.begin(), jac_at_quad.end(), [](vector<double>& local_jac)
  {
    local_jac.resize(1);
    local_jac[0] = 1.;
  });
  return jac_at_quad;
}
*/


} // end namespace Geometry
