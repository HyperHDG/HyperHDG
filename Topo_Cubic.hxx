#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>
#include <array>
#include <vector>

/*!*************************************************************************************************
 * \brief   A namespace containing different classes describing hypergraph topologies.
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
namespace Topology
{  
/*!*************************************************************************************************
 * \brief   Definition of the topology of a hypergraph --- Cubic HyperGraphs.
 *
 * \todo This is not what brief says. It is one special hypergraph.
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
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT       The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template
< unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyEdge_index_t = unsigned int,
  typename hyNode_index_t = hyEdge_index_t >
class Cubic
{
  
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's edges --- Cubic HyperGraph's edges.
   * 
   * \todo    Both private arrays are filled when the hyperedge is constructed. Lazy evaluation
   *          might be an important aspect here. What do you think?
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
    private:
      /*!*******************************************************************************************
       * \brief   Indices of the hypernodes adjacent to the hyperedge.
       *
       * A \c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
       ********************************************************************************************/
      std::array<hyNode_index_t, 2*hyEdge_dimT> hyNode_indices_;
      /*!*******************************************************************************************
       * \brief   Orientation of the hypernode.
       * 
       * \todo    Do we want to change this (cf. detailed description)? This array also does not
       *          have a getter function!
       * 
       * A \c std::array comprising the orientation of each hypernode. In HyperGraph_Cubic, all
       * edges are assumed to have the correct orientation and this array is irrelevant. However, 
       * this is possible to change for different applications.
       ********************************************************************************************/
      std::array<unsigned int, 2*hyEdge_dimT> correct_hyNode_orientation_;
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
      hyEdge(const hyEdge_index_t index, const std::array<unsigned int, space_dimT>& num_elements);
      /*!*******************************************************************************************
       * \brief   Return indices of hypernodes adjacent to the hyperedge.
       *
       * Return a \c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
       *
       * \retval  hypernode_indeices  Topological information on the hyperedge (cf. \c value_type).
       ********************************************************************************************/
      const std::array<hyNode_index_t, 2*hyEdge_dimT>& get_hyNode_indices() const
      { return hyNode_indices_; }
  }; // end of class hyEdge
  
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * \retval  hyEdge_dimT   The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; };
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of the space.
     *
     * \retval  space_dimT       The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; };
  private:
    /*!*********************************************************************************************
     * \brief   Number of elements per spatial dimension.
     *
     * A \c std::array comprising the number of elements in each spatial dimension.
     **********************************************************************************************/
    std::array<unsigned int, space_dimT> num_elements_;
    /*!*********************************************************************************************
     * \brief   Total amount of hyperedges.
     *
     * The number of hyperedges that form the hypergraph. This information is needed to allow to go
     * through all hyperedges and execute some code. The number of hyperedges can be computed from
     * the \c std::array \c num_elements_.
     **********************************************************************************************/
    hyEdge_index_t n_hyEdges_;
    /*!*********************************************************************************************
     * \brief   Total amount of hypernodes.
     *
     * The number of hypernodes that make up the hypergraph. This information is needed to have the
     * appropriate version of a \c HyperNodeFactory. It can be vomputed from the \c std::array
     * \c num_elements_.
     **********************************************************************************************/
    hyNode_index_t n_hyNodes_;
  public:
    /*!*********************************************************************************************
     * \brief   Defines the return value of the class.
     *
     * The \c class \c HyperGraph_Cubic defines the topology of the hypergraph. It "contains" the
     * different hyperedges (that actually are constructed everytime access is needed from e.g. the
     * solver class). Thus, its main purpose is to provide a structure that administrates the
     * hyperedges that are the return value of this structure.
     **********************************************************************************************/
    typedef hyEdge value_type;
    /*!*********************************************************************************************
     * \brief   Defines the value type of input argument for standard constructor.
     *
     * To receive a very general \c AbstractProblem, constructors need to account for the fact that
     * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
     * needs to be constructed. The \c typedef implements the aspect, that a cubic hypergraph
     * topology is by default constructed by a std::vector that contains amounts of elements in the
     * different dimensions.
     **********************************************************************************************/
    typedef std::vector<unsigned int> constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Construct a cubic hypergraph from a \c std::vector.
     *
     * Constructs a hypergraph from a \c std::vector containing the elementens per spatial dimension
     * which is given as input data. If the input vector is shorter that \c space_dimT, the remaining
     * amounts of elemnts are assumed to be equal to zero. If the vector is longer than
     * \c space_dimT, the first \c space_dimT entries are considered only.
     * 
     * \todo    If the vector is too short, an error is thrown in the test program and the behavior
     *          is undefined for Python (most likely an error is thrown, too) at the moment.
     *
     * \param   num_elements    A \c std::vector containing number of elements per dimension.
     **********************************************************************************************/
    Cubic(const constructor_value_type& num_elements);
    /*!*********************************************************************************************
     * \brief   Construct a cubic hypergraph from a \c std::array.
     *
     * \todo    Guido: If possible, this function computes the amount of hyperedges and hypernodes.
     * 
     * Constructs a hypergraph from a \c std::array containing the elementens per spatial dimension
     * which is given as input data. The array has the correct length (as ensured by the involved
     * template parametzer \c space_dimT.
     *
     * \param   num_elements    A \c std::array containing number of elements per spatial dimension.
     **********************************************************************************************/
    Cubic(const std::array<unsigned int, space_dimT>& num_elements);
    /*!*********************************************************************************************
     * \brief   Construct a hypergraph from another hypergraph.
     *
     * \todo    Guido: If possible, this function computes the amount of hyperedges and hypernodes.
     * 
     * Create a (value based) copy of another hypergraph.
     *
     * \param   other           Hypergraph to be copied.
     **********************************************************************************************/
    Cubic(const Cubic<hyEdge_dimT,space_dimT>& other)
    : num_elements_(other.num_elements_), n_hyEdges_(other.n_hyEdges_),
      n_hyNodes_(other.n_hyNodes_) { }
    /*!*********************************************************************************************
     * \brief   Get topological hyperedge of given index.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *
     * This function returns the hyperedge of the given index, i.e., it returns the topological
     * hyperedge (\b not the geometrical information). The topological informatiom comprises the
     * indices of adjacent hypernodes and information about their respective orientations.
     *
     * \param   index           The index of the hyperedge to be returned.
     * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    const value_type operator[](const hyEdge_index_t index) const
    {
      hy_assert( index >= 0 && index < n_hyEdges_ ,
                 "The index of an hyperedge must be non-negative and smaller than the total amount "
                 << "of hyperedges, which is " << n_hyEdges_ << ". Nonetheless, the " << index <<
                 "-th hyperedge is tried to be accessed." );
      return hyEdge(index, num_elements_);
    }
    /*!*********************************************************************************************
     * \brief   Read the array of elements per dimensions.
     *
     * \retval  num_elements    A \c std::array containing the elements in the repective dimension.
     **********************************************************************************************/
    const std::array<unsigned int, space_dimT>& num_elements() const { return num_elements_; }
    /*!*********************************************************************************************
     * \brief   Returns the number of hyperedges making up the hypergraph.
     *
     * \retval  n_hyperedges    The total amount of hyperedges of a hypergraph.
     **********************************************************************************************/
    const hyEdge_index_t n_hyEdges() const { return n_hyEdges_; }
    /*!*********************************************************************************************
     * \brief   Returns the number of hypernodes making up the hypergraph.
     *
     * \retval  n_hypernodes    The total amount of hypernodes of a hypergraph.
     **********************************************************************************************/
    const hyNode_index_t n_hyNodes() const { return n_hyNodes_; }
}; // end of class Cubic


/*
 * HyperEdge functions!
 */


template
< 
  unsigned int space_dimT, typename hyEdge_index_t = unsigned int, 
  typename hyNode_index_t = hyEdge_index_t
>
std::array<hyNode_index_t, 2> line_to_point_index
(const std::array<unsigned int, space_dimT>& num_lines, const hyEdge_index_t index)
{
  hy_assert( num_lines.size() == space_dimT , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  std::array<hyNode_index_t, 2> point_indices;
  point_indices.fill(0);
  
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
    unsigned int helper = 1;
    for (unsigned int dim_n = 0; dim_n < dim_m; ++dim_n)
      helper *= num_lines[dim_n] + 1;
    point_indices[0] += local_indices[dim_m] * helper;
    if (dim_m == orientation)  point_indices[1] += (local_indices[dim_m] + 1) * helper;
    else                       point_indices[1] += local_indices[dim_m] * helper;
  }
  
  return point_indices;
}


template
< 
  unsigned int space_dimT, typename hyEdge_index_t = unsigned int, 
  typename hyNode_index_t = hyEdge_index_t
>
std::array<hyNode_index_t, 4> square_to_line_index
(const std::array<unsigned int, space_dimT>& num_squares, const hyEdge_index_t index)
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
    if (space_dimT == 3)  num_elements_in_direction = num_squares[(dim + orientation) % space_dimT]
                                                       + (dim == 0);
    else                 num_elements_in_direction = num_squares[(dim + orientation) % space_dimT];
    local_indices[(dim + orientation) % space_dimT] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dimT];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No squares should be left any more!" );
  
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
std::array<hyNode_index_t, 6> cube_to_square_index
(const std::array<unsigned int, space_dimT>& num_cubes, const hyEdge_index_t index)
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
  hy_assert( index_helper == 0 , "No cubes should be left any more!" );
  
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


template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyE, typename hyperT>
Cubic<hyEdge_dimT,space_dimT,hyE, hyperT>::hyEdge::
hyEdge(const hyE index, const std::array<unsigned int, space_dimT>& num_elements)
{
  for (unsigned int local_hyNode = 0; local_hyNode < 2 * hyEdge_dimT; ++local_hyNode)
    correct_hyNode_orientation_[local_hyNode] = 1;
  if constexpr ( hyEdge_dimT == 1 )       hyNode_indices_ = line_to_point_index<space_dimT>(num_elements, index);
  else if constexpr ( hyEdge_dimT == 2 )  hyNode_indices_ = square_to_line_index<space_dimT>(num_elements, index);
  else if constexpr ( hyEdge_dimT == 3 )  hyNode_indices_ = cube_to_square_index<space_dimT>(num_elements, index);    
}


/*
 * HyperGraph functions!
 */


template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyE, typename hyperT>
Cubic<hyEdge_dimT,space_dimT,hyE,hyperT>::
Cubic(const std::array<unsigned int, space_dimT>& num_elements)
: num_elements_(num_elements)
{
  static_assert( hyEdge_dimT >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dimT >= hyEdge_dimT, "A domain cannot live within a smaller space!" );
  static_assert( space_dimT <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set n_hyperedges_
  n_hyEdges_ = 1;
  if ( hyEdge_dimT == space_dimT )
    for (unsigned int dim = 0; dim < space_dimT; ++dim)  n_hyEdges_ *= num_elements[dim];
  else if ( hyEdge_dimT == space_dimT - 1 )
  {
    n_hyEdges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyEdges_ += helper;
    }
  }
  else if ( hyEdge_dimT == space_dimT - 2 )
  {
    n_hyEdges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyEdges_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology.");
  hy_assert( n_hyEdges_ > 0 , "An empty hypergraph is being constructed." );
  
  // Set n_hypernodes
  n_hyNodes_ = 1;
  if (hyEdge_dimT == 1)
  {
    n_hyNodes_ *= num_elements[0] + 1;
    if (space_dimT > 1)  n_hyNodes_ *= num_elements[1] + 1;
    if (space_dimT > 2)  n_hyNodes_ *= num_elements[2] + 1;
  }
  else if ( hyEdge_dimT == space_dimT )
  {
    n_hyNodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyNodes_ += helper;
    }
  }
  else if ( hyEdge_dimT == space_dimT - 1 )
  {
    n_hyNodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyNodes_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hyNodes_ > 0 , "An empty hypergraph is being constructed." );
}

template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename hyE, typename hyT>
Cubic<hyEdge_dimT,space_dimT,hyE,hyT>::
Cubic(const constructor_value_type& num_elements)
{
  for (unsigned int dim = 0; dim < space_dimT; ++dim) num_elements_[dim] = num_elements[dim];
  
  static_assert( hyEdge_dimT >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dimT >= hyEdge_dimT, "A domain cannot live within a smaller space!" );
  static_assert( space_dimT <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set n_hyperedges_
  n_hyEdges_ = 1;
  if ( hyEdge_dimT == space_dimT )
    for (unsigned int dim = 0; dim < space_dimT; ++dim)  n_hyEdges_ *= num_elements[dim];
  else if ( hyEdge_dimT == space_dimT - 1 )
  {
    n_hyEdges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyEdges_ += helper;
    }
  }
  else if ( hyEdge_dimT == space_dimT - 2 )
  {
    n_hyEdges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyEdges_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hyEdges_ > 0 , "An empty hypergraph is being constructed." );
  
  // Set n_hypernodes
  n_hyNodes_ = 1;
  if (hyEdge_dimT == 1)
  {
    n_hyNodes_ *= num_elements[0] + 1;
    if (space_dimT > 1)  n_hyNodes_ *= num_elements[1] + 1;
    if (space_dimT > 2)  n_hyNodes_ *= num_elements[2] + 1;
  }
  else if ( hyEdge_dimT == space_dimT )
  {
    n_hyNodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyNodes_ += helper;
    }
  }
  else if ( hyEdge_dimT == space_dimT - 1 )
  {
    n_hyNodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dimT; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dimT; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyNodes_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hyNodes_ > 0 , "An empty hypergraph is being constructed." );
}







} // end of namespace Topology
