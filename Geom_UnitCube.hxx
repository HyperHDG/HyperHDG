#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/Hypercube.hxx>
#include <HyperHDG/DenseLA.hxx>
#include <HyperHDG/Topology/Cubic.hxx>
#include <HyperHDG/Mapping/Linear.hxx>

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
 * \todo    CONSIDER WRAPPER CONSTRUCTOR!
 * 
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
template 
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, typename pt_coord_t = double,
  template < unsigned int, unsigned int, typename >  typename mapping_tM = Mapping::Linear,
  unsigned int hyEdge_dimTM = hyEdge_dimT, unsigned int space_dimTM = space_dimT,
  typename pt_coord_tM = pt_coord_t, typename hyEdge_index_t = unsigned int
>
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
    /*!*********************************************************************************************
     * \brief   The mapping type is \c mapping_tt with given template parameters.
     **********************************************************************************************/
    using mapping_t = mapping_tM<hyEdge_dimTM,space_dimTM,pt_coord_tM>;
    private:
      /*!*******************************************************************************************
       * \brief   Points adjacent to the hyperedge.
       *
       * \todo    In the long run, do not store these.
       *
       * An array comprising the vertices (points) of a cubic hyperedge.
       ********************************************************************************************/
      std::array<Point<space_dimT>, Hypercube<hyEdge_dimT>::n_vertices()> points_;
      /*!*******************************************************************************************
       * \brief   Fill array of vertices of hyEdge.
       ********************************************************************************************/
      template<unsigned int hyEdge_dimTT, unsigned int space_dimTT>
      unsigned int fill_points
      ( 
        unsigned int index, 
        const tpcc_elem_t<hyEdge_dimTT,space_dimTT>& elem, const UnitCube& geometry
      )
      {
        if constexpr ( hyEdge_dimTT == 0 )
        {
          Point<space_dimT> pt;
          for (unsigned int dim = 0; dim < space_dimTT; ++dim)
            pt[dim] = exterior_coordinate<hyEdge_dimTT, space_dimT>(elem, dim)
                        / geometry.num_elements_[dim];
          points_[index++] = pt;
        }
        else
          for (unsigned int i = 0; i < 2; ++i)
            index = fill_points<hyEdge_dimTT-1,space_dimTT>
                      ( index, get_face<hyEdge_dimTT, space_dimTT>(elem, i), geometry );
        return index;
      }
      /*!*******************************************************************************************
       * \brief   Hold an instance of a mapping type to be able to calculate normals and so on.
       ********************************************************************************************/
      std::shared_ptr<mapping_t> mapping;
      /*!*******************************************************************************************
       * \brief   Prepare lazy evaluation of mapping.
       ********************************************************************************************/
      inline void generate_mapping_if_needed()
      {
        if (mapping) return;
        Point<space_dimT,pt_coord_t> translation = (Point<space_dimT,pt_coord_t>) point(0);
        SmallMat<space_dimT,hyEdge_dimT,pt_coord_t> matrix;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          matrix.set_column(dim, (Point<space_dimT,pt_coord_t>) point(1<<dim) - translation);
        mapping = std::make_shared<mapping_t>(translation, matrix);
      }
      /*!*******************************************************************************************
       * \brief   Hold an instance of a matrix containing local normals.
       ********************************************************************************************/
      std::shared_ptr< SmallSquareMat<hyEdge_dimT,pt_coord_t> > local_normals_;
      /*!*******************************************************************************************
       * \brief   Hold an instance of a matrix containing inner and outer normals.
       ********************************************************************************************/
      std::shared_ptr< SmallSquareMat<space_dimT,pt_coord_t> > global_normals_;
    public:
      /*!*******************************************************************************************
       * \brief   Returns dimension of the hyperedge.
       ********************************************************************************************/
      static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
      /*!*******************************************************************************************
       * \brief   Returns dimension of the surrounding space.
       ********************************************************************************************/
      static constexpr unsigned int space_dim() { return space_dimT; }
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
        tpcc_elem_t<hyEdge_dimT, space_dimT> elem 
          = get_element<hyEdge_dimT, space_dimT, hyEdge_index_t>(geometry.tpcc_elements_, index);
        fill_points<hyEdge_dimT,space_dimT>(0, elem, geometry);
      }
      /*!*******************************************************************************************
       * \brief   Return vertex of specified index of a hyperedge.
       *
       * Return a \c Point describing the position of a vertex of a hyperedge.
       *
       * \retval  point           Point/Vertex of the hyperedge.
       ********************************************************************************************/
      Point<space_dimT> point(const unsigned int index) const  { return points_[index]; }


      /*!*******************************************************************************************
       * \brief   Map n_vec points from reference to physical element.
       ********************************************************************************************/
      template < unsigned int n_vec >
      SmallMat<space_dimT, n_vec, pt_coord_t> map_ref_to_phys
      (const SmallMat<hyEdge_dimT, n_vec, pt_coord_t>& pts)
      {
        for (unsigned int i = 0; i < pts.size(); ++i)
          hy_assert( pts[i] >= 0. && pts[i] <= 1. ,
                     "Point must lie in reference square!");
        generate_mapping_if_needed();
        return mapping->map_reference_to_physical(pts);
      }
      /*!*******************************************************************************************
       * \brief   Return matrix column of the affine-linear transformation.
       * 
       * \param   index   Index of the matrix column to be returned.
       * \retval  column  The specified matrix column.
       ********************************************************************************************/
      SmallVec<space_dimT, pt_coord_t> span_vec(const unsigned int index)
      {
        hy_assert( index < hyEdge_dimT,
                   "There are only " << hyEdge_dimT << " spanning vectors." );
        generate_mapping_if_needed();
        return mapping->get_column(index);
      }
      /*!*******************************************************************************************
       * \brief   Return reduced matrix R of the QR decomposition.
       ********************************************************************************************/
      const SmallSquareMat<hyEdge_dimT,pt_coord_t>& mat_r()
      {
        generate_mapping_if_needed();
        return mapping->mat_r();
      }
      /*!*******************************************************************************************
       * \brief   Return Haussdorff/Lebesque measure of the hyperedge.
       ********************************************************************************************/
      pt_coord_t area()
      {
        generate_mapping_if_needed();
        return std::abs(mapping->functional_determinant_hyEdge());
      }
      /*!*******************************************************************************************
       * \brief   Return Haussdorff measure of the specified hypernode.
       ********************************************************************************************/
      pt_coord_t face_area(const unsigned int  index)
      {
        hy_assert( index < 2 * hyEdge_dimT ,
                   "A hyperedge has 2 * dim(hyEdge) faces." );
        generate_mapping_if_needed();
        return std::abs(mapping->functional_determinant_hyNode(index / 2));
      }
      /*!*******************************************************************************************
       * \brief   Return local normal of given index.
       *
       * Return outer unit normal with respect to the hypernode which is spanned by the vectors
       * spanning the phyiscal element, orthogonally projected to a hyEdge_dimT dimensiona space,
       * but the vector of the given index. This is an element of the same dimension as the 
       * reference element.
       ********************************************************************************************/
      Point<hyEdge_dimT,pt_coord_t> local_normal(const unsigned int index)
      {
        hy_assert( index < 2 * hyEdge_dimT ,
                   "A hyperedge has 2 * dim(hyEdge) inner normals." );
        generate_mapping_if_needed();
        if (!local_normals_)
          local_normals_ = std::make_shared< SmallSquareMat<hyEdge_dimT,pt_coord_t> > ();
        Point<hyEdge_dimT,pt_coord_t> normal = local_normals_->get_column(index / 2);
        if (norm_2(normal) < 0.5)
        {
          normal = mapping->local_normal(index / 2);
          local_normals_->set_column(index / 2, normal);
        }
        if (index % 2 == 1)  normal *= -1.;
        return normal;
      }
      /*!*******************************************************************************************
       * \brief   Return inner normal of given index.
       *
       * Return outer unit normal with respect to the hypernode which is spanned by all vectors
       * spanning the phyiscal element, but the vector of the given index. The vector has to be in 
       * the span of the columns of the local transformation matrix. This is an element of the same
       * dimension as the full space.
       ********************************************************************************************/
      Point<space_dimT,pt_coord_t> inner_normal(const unsigned int index)
      {
        hy_assert( index < 2 * hyEdge_dimT ,
                   "A hyperedge has 2 * dim(hyEdge) inner normals." );
        generate_mapping_if_needed();
        if (!global_normals_)
          global_normals_ = std::make_shared< SmallSquareMat<space_dimT,pt_coord_t> > ();
        Point<space_dimT,pt_coord_t> normal = global_normals_->get_column(index / 2);
        if (norm_2(normal) < 0.5)
        {
          normal = mapping->inner_normal(index / 2);
          global_normals_->set_column(index / 2, normal);
        }
        if (index % 2 == 1)  normal *= -1.;
        return normal;
      }
      /*!*******************************************************************************************
       * \brief   Return outer normal of given index.
       *
       * Return unit normal with respect to the hyperedge within the full space.
       ********************************************************************************************/
      Point<space_dimT,pt_coord_t> outer_normal(const unsigned int index)
      {
        hy_assert( index < space_dimT - hyEdge_dimT ,
                   "This function returns one of the dim(space) - dim(hyEdge) orthonormal vectors "
                   << "which are orthogonal to the hyperedge." );
        generate_mapping_if_needed();
        if (!global_normals_)
          global_normals_ = std::make_shared< SmallSquareMat<space_dimT,pt_coord_t> > ();
        Point<space_dimT,pt_coord_t> normal = global_normals_->get_column(hyEdge_dimT + index);
        if (norm_2(normal) < 0.5)
        {
          normal = mapping->outer_normal(index);
          global_normals_->set_column(hyEdge_dimT + index, normal);
        }
        return normal;
      }
      /*!*******************************************************************************************
       * \brief   Return lexicographically ordered equidistant tensorial point of given index.
       ********************************************************************************************/
      template<unsigned int n_sub_points, typename one_dim_float_t>
      Point<space_dimT,pt_coord_t> lexicographic
      ( unsigned int index, const std::array<one_dim_float_t, n_sub_points>& points_1d )
      {
        static_assert( n_sub_points > 0 , "No subpoints do not make sense!" );
        hy_assert( index < std::pow(n_sub_points, hyEdge_dimT) ,
                   "The index must niot exceed the number of prescribed lexicographic points." );
        generate_mapping_if_needed();
        Point<hyEdge_dimT,pt_coord_t> pt;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          pt[dim] = (pt_coord_t) points_1d[index % n_sub_points]; 
          index /= n_sub_points;
        }
        return mapping->map_reference_to_physical(pt);
      }
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
    tpcc_t<hyEdge_dimT, space_dimT, hyEdge_index_t> tpcc_elements_;
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
    : tpcc_elements_(num_elements_)
    { 
      for (unsigned int dim = 0; dim < space_dimT; ++dim) num_elements_[dim] = num_elements[dim];
      tpcc_elements_ = create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements_);
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
      tpcc_elements_(create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(other.num_elements()))
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

} // end namespace Geometry
