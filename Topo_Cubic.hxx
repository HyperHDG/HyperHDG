#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/Wrapper/TPCC.hxx>
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
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
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
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
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
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
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
//      std::array<unsigned int, 2*hyEdge_dimT> correct_hyNode_orientation_;
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
      hyEdge(const hyEdge_index_t index, const Cubic& topology)
      {
        tpcc_elem_t<hyEdge_dimT, space_dimT> elem 
          = get_element<hyEdge_dimT, space_dimT,unsigned int>(topology.tpcc_elements_, index);
        for (unsigned int i = 0; i < hyNode_indices_.size(); ++i)
        {
          tpcc_elem_t<hyEdge_dimT-1, space_dimT> face = get_face<hyEdge_dimT, space_dimT>(elem, i);
          hyNode_indices_[i]
            = get_index<hyEdge_dimT-1, space_dimT,unsigned int>(topology.tpcc_faces_, face);
        }
      }
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
     * \brief   Tensor product chain complex.
     **********************************************************************************************/
    tpcc_t<hyEdge_dimT, space_dimT, hyNode_index_t> tpcc_elements_;
    tpcc_t<hyEdge_dimT-1, space_dimT, hyNode_index_t> tpcc_faces_;
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
    Cubic(const constructor_value_type& num_elements)
    : tpcc_elements_(num_elements_), tpcc_faces_(num_elements_)
    {
      hy_assert( num_elements.size() == space_dimT ,
                 "Incorrect size of number of elements specifications!" );
      for (unsigned int dim = 0; dim < space_dimT; ++dim)  num_elements_[dim] = num_elements[dim];
      tpcc_elements_ = create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements_);
      tpcc_faces_ = tpcc_faces< hyEdge_dimT, space_dimT, hyEdge_index_t >(tpcc_elements_);
      n_hyEdges_ = n_elements< hyEdge_dimT, space_dimT, hyEdge_index_t >(tpcc_elements_);
      n_hyNodes_ = n_elements< hyEdge_dimT-1, space_dimT, hyEdge_index_t >(tpcc_faces_);
    }
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
    Cubic(const std::array<unsigned int, space_dimT>& num_elements)
    : num_elements_(num_elements),
      tpcc_elements_(create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements)),
      tpcc_faces_(tpcc_faces< hyEdge_dimT, space_dimT, hyEdge_index_t >(tpcc_elements_)),
      n_hyEdges_(n_elements< hyEdge_dimT, space_dimT, hyEdge_index_t >(tpcc_elements_)), n_hyNodes_(n_elements< hyEdge_dimT-1, space_dimT, hyEdge_index_t >(tpcc_faces_))
    { }
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
    : num_elements_(other.num_elements_), tpcc_elements_(other.tpcc_elements_),
      tpcc_faces_(other.tpcc_faces_), n_hyEdges_(other.n_hyEdges_), n_hyNodes_(other.n_hyNodes_) { }
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
      return hyEdge(index, *this);
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

} // end of namespace Topology
