#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/mapping/linear.hxx>
#include <HyperHDG/topology/file.hxx>

namespace Geometry
{
/*!*************************************************************************************************
 * \brief   Hypergraph geometry based on an input file.
 *
 * The geometry class File is a set of hyperedges. Each of these tensorial hyperedges is represented
 * by its vertices (given within the file). For consistency, it is assumed that the vertices are
 * given in lexicographical order. Also the hypernodes are assumed to be given in lexicographical
 * order to ensure that geometry and topology of all hyperedges fit.
 *
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  vectorT         The typename of the large vector type. Defaults to std::vector.
 * \tparam  pointT          The typename of a Point class. Defaults to \c Point<space_dimT, float>.
 * \tparam  mapping_tM      The mapping which should depened on the aforementioned three template
 *                          arguments. Default is linear mapping.
 * \tparam  hyEdge_dimTM    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 *                          Template parameter for the mapping which defaults to hyEdge_dimT.
 * \tparam  space_dimTM     The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 *                          Template parameter for the mapping which defaults to space_dimT.
 * \tparam  pt_coord_tM     The floating point type in which points/vertices are represented.
 *                          Default is the floating point of the point class.
 *                          Template parameter for the mapping which defaults to pt_coord_t.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is \c hyNode_index_t.
 * \tparam  pt_index_t      The index type of points. Default is \c hyNode_index_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          template <typename...> typename vectorT = std::vector,
          typename pointT = Point<space_dimT, double>,
          template <unsigned int, unsigned int, typename> typename mapping_tM = Mapping::Linear,
          unsigned int hyEdge_dimTM = hyEdge_dimT,
          unsigned int space_dimTM = space_dimT,
          typename pt_coord_tM = typename pointT::value_type,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
class File
{
  /*!***********************************************************************************************
   * \brief   The point type is the floating point type of \c pointT.
   ************************************************************************************************/
  using pt_coord_t = typename pointT::value_type;
  /*!***********************************************************************************************
   * \brief   The mapping type is \c mapping_tt with given template parameters.
   ************************************************************************************************/
  using mapping_t = mapping_tM<hyEdge_dimTM, space_dimTM, pt_coord_tM>;
  /*!***********************************************************************************************
   * \brief   Definition of the geometry of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
   public:
    /*!*********************************************************************************************
     * \brief   Return dimension of the hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Return dimension of the surrounding space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; }

   private:
    /*!*********************************************************************************************
     * \brief   Reference to parent hypergraph.
     **********************************************************************************************/
    const File& hyGraph_geometry_;
    /*!*********************************************************************************************
     * \brief   Index of the hyperedge within the hypergraph.
     **********************************************************************************************/
    const hyEdge_index_t index_;
    /*!*********************************************************************************************
     * \brief   Hold an instance of a mapping type to be able to calculate normals and so on.
     **********************************************************************************************/
    mapping_t mapping;
    /*!*********************************************************************************************
     * \brief   Return vertex of specified index of a hyperedge.
     **********************************************************************************************/
    pointT point(const unsigned int pt_index) const
    {
      return hyGraph_geometry_.domain_info_
        .points[hyGraph_geometry_.domain_info_
                  .points_hyEdge[index_ / hyGraph_geometry_.n_loc_ref_elem][pt_index]];
    }

   public:
    /*!*********************************************************************************************
     * \brief   Construct hyperedge from hypergraph and index.
     **********************************************************************************************/
    hyEdge(const File& hyGraph_geometry, const hyEdge_index_t index)
    : hyGraph_geometry_(hyGraph_geometry), index_(index)
    {
      Point<space_dimT, pt_coord_t> translation = point(0);
      SmallMat<space_dimT, hyEdge_dimT, pt_coord_t> matrix;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        matrix.set_column(dim, point(1 << dim) - translation);

      Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT> elem = Wrapper::get_element(
        hyGraph_geometry.tpcc_ref_elem_, index % hyGraph_geometry_.n_loc_ref_elem);
      matrix /= (pt_coord_t)hyGraph_geometry_.n_subintervals_;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        translation += (pt_coord_t)Wrapper::interior_coordinate(elem, dim) * matrix.get_column(dim);

      mapping = mapping_t(translation, matrix);
    }
    /*!*********************************************************************************************
     * \brief   Map n_vec points from reference to physical element.
     *
     * \param   points        Matrix whose columns consist of the points to be mapped.
     * \retval  phy_points    Matrix whose columns consist of the mapped points.
     **********************************************************************************************/
    template <unsigned int n_vec>
    SmallMat<space_dimT, n_vec, pt_coord_t> map_ref_to_phys(
      const SmallMat<hyEdge_dimT, n_vec, pt_coord_t>& points) const
    {
      for (unsigned int i = 0; i < points.size(); ++i)
        hy_assert(points[i] >= 0. && points[i] <= 1., "Point must lie in reference square!");
      return mapping.map_reference_to_physical(points);
    }
    /*!*********************************************************************************************
     * \brief   Map n_vec points from reference to physical element.
     *
     * \param   points        Matrix whose columns consist of the points to be mapped.
     * \retval  points        Matrix whose columns consist of the mapped points.
     **********************************************************************************************/
    template <unsigned int n_vec>
    SmallMat<space_dimT, n_vec, pt_coord_t>& map_ref_to_phys(
      SmallMat<space_dimT, n_vec, pt_coord_t>& points) const
    {
      for (unsigned int i = 0; i < points.size(); ++i)
        hy_assert(points[i] >= 0. && points[i] <= 1., "Point must lie in reference square!");
      points = mapping.map_reference_to_physical(points);
      return points;
    }
    /*!*********************************************************************************************
     * \brief   Return matrix column of the affine-linear transformation.
     *
     * \param   index   Index of the matrix column to be returned.
     * \retval  column  The specified matrix column.
     **********************************************************************************************/
    SmallVec<space_dimT, pt_coord_t> span_vec(const unsigned int index) const
    {
      hy_assert(index < hyEdge_dimT, "There are only " << hyEdge_dimT << " spanning vectors.");
      return mapping.get_column(index);
    }
    /*!*********************************************************************************************
     * \brief   Return reduced matrix R of the QR decomposition.
     **********************************************************************************************/
    const SmallSquareMat<hyEdge_dimT, pt_coord_t>& mat_r() const { return mapping.mat_r(); }
    /*!*********************************************************************************************
     * \brief   Return matrix Q of the QR decomposition of the linear transoformation.
     **********************************************************************************************/
    const SmallSquareMat<space_dimT, pt_coord_t> mat_q() const { return mapping.mat_q(); }
    /*!*********************************************************************************************
     * \brief   Return Haussdorff/Lebesque measure of the hyperedge.
     **********************************************************************************************/
    pt_coord_t area() const { return std::abs(mapping.functional_determinant_hyEdge()); }
    /*!*********************************************************************************************
     * \brief   Return Haussdorff measure of the specified hypernode.
     **********************************************************************************************/
    pt_coord_t face_area(const unsigned int index) const
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) faces.");
      return std::abs(mapping.functional_determinant_hyNode(index / 2));
    }
    /*!*********************************************************************************************
     * \brief   Return local normal of given index.
     *
     * Return outer unit normal with respect to the hypernode which is spanned by the vectors
     * spanning the phyiscal element, orthogonally projected to a hyEdge_dimT dimensional space,
     * but the vector of the given index. This is an element of the same dimension as the
     * reference element.
     **********************************************************************************************/
    Point<hyEdge_dimT, pt_coord_t> local_normal(const unsigned int index) const
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) inner normals.");
      Point<hyEdge_dimT, pt_coord_t> normal = mapping.local_normal(index / 2);
      if (index % 2 == 1)
        normal *= -1.;
      return normal;
    }
    /*!*********************************************************************************************
     * \brief   Return inner normal of given index.
     *
     * Return outer unit normal with respect to the hypernode which is spanned by all vectors
     * spanning the phyiscal element, but the vector of the given index. The vector has to be in
     * the span of the columns of the local transformation matrix. This is an element of the same
     * dimension as the full space.
     **********************************************************************************************/
    Point<space_dimT, pt_coord_t> inner_normal(const unsigned int index) const
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) inner normals.");
      Point<space_dimT, pt_coord_t> normal = mapping.inner_normal(index / 2);
      if (index % 2 == 1)
        normal *= -1.;
      return normal;
    }
    /*!*********************************************************************************************
     * \brief   Return barycenter of face of given index.
     **********************************************************************************************/
    Point<space_dimT, pt_coord_t> face_barycenter(const unsigned int index) const
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) inner normals.");
      Point<hyEdge_dimT, pt_coord_t> local_center(0.5);
      local_center[index / 2] = index % 2 ? 1. : 0.;
      return map_ref_to_phys(local_center);
    }
    /*!*********************************************************************************************
     * \brief   Return outer normal of given index.
     *
     * Return unit normal with respect to the hyperedge within the full space.
     **********************************************************************************************/
    Point<space_dimT, pt_coord_t> outer_normal(const unsigned int index) const
    {
      hy_assert(index < space_dimT - hyEdge_dimT,
                "This function returns one of the dim(space) - dim(hyEdge) orthonormal vectors "
                  << "which are orthogonal to the hyperedge.");
      return mapping.outer_normal(index);
    }
    /*!*********************************************************************************************
     * \brief   Return lexicographically ordered equidistant tensorial point of given index.
     **********************************************************************************************/
    template <unsigned int n_sub_points, typename one_dim_float_t>
    Point<space_dimT, pt_coord_t> lexicographic(
      unsigned int index,
      const SmallVec<n_sub_points, one_dim_float_t>& points_1d) const
    {
      static_assert(n_sub_points > 0, "No subpoints do not make sense!");
      hy_assert(index < std::pow(n_sub_points, hyEdge_dimT),
                "The index must not exceed the number of prescribed lexicographic points.");
      Point<hyEdge_dimT, pt_coord_t> pt;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        pt[dim] = (pt_coord_t)points_1d[index % n_sub_points];
        index /= n_sub_points;
      }
      return mapping.map_reference_to_physical(pt);
    }
    /*!*********************************************************************************************
     * \brief   Return equidistant tensorial point of given index on a given boundary (slightly
     *          moved away from the boundary using boundary_scale), ordered lexicographically.
     **********************************************************************************************/
    template <unsigned int n_sub_points, typename one_dim_float_t>
    Point<space_dimT, pt_coord_t> boundary_lexicographic(
      unsigned int index,
      unsigned int boundary_number,
      float boundary_scale,
      const SmallVec<n_sub_points, one_dim_float_t>& points_1d) const
    {
      static_assert(n_sub_points > 0, "No subpoints do not make sense!");
      hy_assert(index < std::pow(n_sub_points, hyEdge_dimT - 1) * hyEdge_dimT * 2,
                "The index must not exceed the number of prescribed lexicographic points.");
      Point<hyEdge_dimT - 1, pt_coord_t> subpt;
      index = index - std::pow(n_sub_points, hyEdge_dimT - 1) * boundary_number;
      for (unsigned int subdim = 0; subdim < hyEdge_dimT - 1; ++subdim)
      {
        subpt[subdim] = (pt_coord_t)points_1d[index % n_sub_points];
        index /= n_sub_points;
      }
      Point<hyEdge_dimT, pt_coord_t> pt;
      unsigned int subd = 0;
      for (unsigned int d = 0; d < hyEdge_dimT; d++)
      {
        if (boundary_number / 2 == d)
          pt[d] = boundary_scale * (boundary_number % 2 - 0.5) + 0.5;
        else
          pt[d] = subpt[subd++];
      }
      return mapping.map_reference_to_physical(pt);
    }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT   The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT       The dimension of the space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; }

 private:
  /*!***********************************************************************************************
   * \brief   Domain Info containing all the information of the hypergraph (cf. ReadDomain.hxx).
   ************************************************************************************************/
  const DomainInfo<hyEdge_dimT,
                   space_dimT,
                   vectorT,
                   pointT,
                   hyEdge_index_t,
                   hyNode_index_t,
                   pt_index_t>& domain_info_;

  /*!***********************************************************************************************
   * \brief   Refinment level corresponds to number of subintervals per dimension.
   ************************************************************************************************/
  unsigned int n_subintervals_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for refined elements.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t> tpcc_ref_elem_;
  /*!***********************************************************************************************
   * \brief   Number of refined elements per element.
   ************************************************************************************************/
  unsigned int n_loc_ref_elem;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class.
   *
   * The \c class \c Geometry::File defines the geometry of the hypergraph. It "contains" the
   * different hyperedges (that actually are constructed everytime access is needed from e.g. the
   * solver class). Thus, its main purpose is to provide a structure that administrates the
   * hyperedges that are the return value of this structure.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Define the value type of input argument for standard constructor.
   ************************************************************************************************/
  typedef Topology::
    File<hyEdge_dimT, space_dimT, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
      constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a hypergraph which is read from a file.
   *
   * Constructs a hypergraph from a \c Topology::File containing the elementens per spatial
   * dimension which is given as by its topology.
   *
   * \param   topology     The topology of the hypergraph that has the geometry of the unit cube.
   ************************************************************************************************/
  File(const constructor_value_type& topology)
  : domain_info_(topology.domain_info()),
    n_subintervals_(topology.get_refinement()),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervals_))),
    n_loc_ref_elem(Hypercube<hyEdge_dimT>::pow(n_subintervals_))
  {
  }
  /*!***********************************************************************************************
   * \brief   Get geometrical hyperedge of given index.
   *
   * This function returns the hyperedge of the given index, i.e., it returns the geometrical
   * hyperedge (\b not the topological information). The geometrical informatiom comprises the
   * indices of adjacent vertices (i.e. points) and information about their respective positions.
   * It is equivalent to \c get_hyEdge.
   *
   * \param   index       The index of the hyperedge to be returned.
   * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  value_type operator[](const hyEdge_index_t index) const { return get_hyEdge(index); }
  /*!***********************************************************************************************
   * \brief   Get geometrical hyperedge of given index.
   *
   * This function returns the hyperedge of the given index, i.e., it returns the geometrical
   * hyperedge (\b not the topological information). The geometrical informatiom comprises the
   * indices of adjacent vertices (i.e. points) and information about their respective positions.
   * It is equivalent to \c operator[].
   *
   * \param   index       The index of the hyperedge to be returned.
   * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  value_type get_hyEdge(const hyEdge_index_t index) const
  {
    hy_assert(index < domain_info_.n_hyEdges * Wrapper::n_elements(tpcc_ref_elem_) && index >= 0,
              "Index must be non-negative and smaller than "
                << domain_info_.n_hyEdges << " (which is the amount of hyperedges). It was "
                << index << "!");
    return hyEdge(*this, index);
  }
  /*!***********************************************************************************************
   * \brief   Return the refinement level (equal to number of subintervals).
   ************************************************************************************************/
  unsigned int get_refinement() const { return n_subintervals_; }
  /*!***********************************************************************************************
   * \brief   Set the refinement level (equal to number of subintervals).
   ************************************************************************************************/
  void set_refinement(unsigned int level)
  {
    n_subintervals_ = level;
    tpcc_ref_elem_ =
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervals_));
    n_loc_ref_elem = Hypercube<hyEdge_dimT>::pow(n_subintervals_);
  }
};  // end class File

}  // end namespace Geometry
