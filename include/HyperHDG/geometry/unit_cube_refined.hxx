#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/cubic_refined.hxx>

#include <algorithm>

namespace Geometry
{
/*!*************************************************************************************************
 * \brief   Define the geometry of a unit cube that is subdivided into several orthotopes.
 *
 * This class defines the geometry of a unit cube which is subdivided into several orthotopes. The
 * number of orthotopes in each spatial dimension is defined by the vector that serves as the
 * argument of the constructor.
 *
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes, ....
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  pt_coord_t      The floating point type in which the coordinates of vertices are given.
 *                          Defaults to double.
 * \tparam  ConstructorVecT The vector/array type of the constructor.
 *                          Defaults to SmallVec<space_dimT, unsigned int>.
 * \tparam  hyEdge_index_t  The integer type in which the indices of hyperedges are represented.
 *                          Defaults to unsigne int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          unsigned int n_subintervalsT = 1,
          typename pt_coord_t = double,
          typename ConstructorVecT = SmallVec<space_dimT, unsigned int>,
          typename hyEdge_index_t = unsigned int>
class UnitCubeRefined
{
  /*!***********************************************************************************************
   * \brief   Definition of the geometry of a single hyperedge.
   *
   * A hyperege is uniquely defined by a translation, the spanning dimensions, and the lengths of
   * the lines along the these spanning dimensions.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
   private:
    /*!*********************************************************************************************
     * \brief   Translation vector of the origin.
     **********************************************************************************************/
    Point<space_dimT, pt_coord_t> translation;
    /*!*********************************************************************************************
     * \brief   The dimensions in wihch the orthotope evolves.
     **********************************************************************************************/
    SmallVec<hyEdge_dimT, unsigned int> dim_indices;
    /*!*********************************************************************************************
     * \brief   Lengths of the orthotope's lines which are parallel to the \c dim_indices axes.
     **********************************************************************************************/
    SmallVec<hyEdge_dimT, pt_coord_t> char_length;

    /*!*********************************************************************************************
     * \brief   Fill data of hyEdge.
     *
     * Recursive function that fills \c translation, \c dim_indices, and \¢ char_length. To do so,
     * it considers a single hyperedge and takes its zeroth and first faces, and the zeroth and
     * first faces of these faces, and so on to determine the positions of the translation and the
     * other spanning points of the interface.
     *
     * \tparam  hyEdge_dimTT  Dimension of the hyperedge or face that is considered.
     * \param   index         Variable that indicates which point is investigated.
     *                        Supposed to be 0, initially.
     * \param   elem          Hyperedge or face which is considered.
     * \param   geometry      The overall unit cube of which the initial elem is a hyperedge of.
     * \retval  index         Index that is used for further (recursive) executions of the function.
     **********************************************************************************************/
    template <unsigned int hyEdge_dimTT>
    unsigned int fill_data(unsigned int index,
                           const Wrapper::tpcc_elem_t<hyEdge_dimTT, space_dimT>& elem,
                           const UnitCubeRefined& geometry)
    {
      if constexpr (hyEdge_dimTT == 0)
      {
        if (index == 0)
          for (unsigned int dim_pt = 0; dim_pt < space_dimT; ++dim_pt)
          {
            unsigned int ext_dim = Wrapper::exterior_direction(elem, dim_pt);
            translation[ext_dim] = (pt_coord_t)Wrapper::exterior_coordinate(elem, dim_pt) /
                                   (pt_coord_t)geometry.n_elements_[ext_dim];
            hy_assert(0. <= translation[ext_dim] && translation[ext_dim] <= 1.,
                      "The unit cube has only these cooridnates.");
          }
        unsigned int dim = 0;
        for (; dim < hyEdge_dimT && char_length[dim] != 0.; ++dim)
          ;
        if (index == (unsigned int)1 << dim && char_length[dim] == 0.)
          for (unsigned int dim_pt = 0; dim_pt < space_dimT; ++dim_pt)
          {
            unsigned int ext_dim = Wrapper::exterior_direction(elem, dim_pt);
            pt_coord_t helper = (pt_coord_t)Wrapper::exterior_coordinate(elem, dim_pt) /
                                (pt_coord_t)geometry.n_elements_[ext_dim];
            hy_assert(0. <= helper && helper <= 1., "The unit cube has only these cooridnates.");
            if (helper != translation[ext_dim])
            {
              char_length[dim] = helper - translation[ext_dim];
              dim_indices[dim] = ext_dim;
              break;
            }
          }
        ++index;
      }
      else
      {
        for (unsigned int i = 2 * hyEdge_dimTT - 2; i < 2 * hyEdge_dimTT; ++i)
          index = fill_data<hyEdge_dimTT - 1>(index, Wrapper::get_face(elem, i), geometry);
      }
      return index;
    }

    void adapt_data(const Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT>& elem)
    {
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
      {
        char_length[dim] /= (pt_coord_t)n_subintervalsT;
        translation[dim_indices[dim]] +=
          (pt_coord_t)Wrapper::interior_coordinate(elem, dim) * char_length[dim];
      }
    }

   public:
    /*!*********************************************************************************************
     * \brief   Return dimension of the hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Return dimension of the surrounding space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; }
    /*!*********************************************************************************************
     * \brief   Construct a orthotopic hyperedge from its index and the surrounding unit hypercube.
     *
     * \param   index         The index of the hyperedge to be created.
     * \param   geometry      The surrounding unit hypercube.
     **********************************************************************************************/
    hyEdge(const hyEdge_index_t index, const UnitCubeRefined& geometry)
    {
      Wrapper::tpcc_elem_t<hyEdge_dimT, space_dimT> elem =
        Wrapper::get_element(geometry.tpcc_elements_, index / geometry.n_loc_ref_elem);
      Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT> loc_elem =
        Wrapper::get_element(geometry.tpcc_ref_elem_, index % geometry.n_loc_ref_elem);
      fill_data<hyEdge_dimT>(0, elem, geometry);
      adapt_data(loc_elem);
    }

    /*!*********************************************************************************************
     * \brief   Map n_vec points from reference to physical element.
     *
     * \param   points        Matrix whose columns consist of the points to be mapped.
     * \retval  phy_points    Matrix whose columns consist of the mapped points.
     **********************************************************************************************/
    template <unsigned int n_vec>
    SmallMat<space_dimT, n_vec, pt_coord_t> map_ref_to_phys(
      const SmallMat<hyEdge_dimT, n_vec, pt_coord_t>& points)
    {
      for (unsigned int i = 0; i < points.size(); ++i)
        hy_assert(points[i] >= 0. && points[i] <= 1., "Point must lie in reference square!");

      SmallMat<space_dimT, n_vec, pt_coord_t> phy_points =
        rep_mat<space_dimT, n_vec, pt_coord_t>(translation);
      for (unsigned int j = 0; j < n_vec; ++j)
        for (unsigned int i = 0; i < hyEdge_dimT; ++i)
          phy_points(dim_indices[i], j) += points(i, j) * char_length[i];

      return phy_points;
    }
    /*!*********************************************************************************************
     * \brief   Map n_vec points from reference to physical element.
     *
     * \param   points        Matrix whose columns consist of the points to be mapped.
     * \retval  points        Matrix whose columns consist of the mapped points.
     **********************************************************************************************/
    template <unsigned int n_vec>
    SmallMat<space_dimT, n_vec, pt_coord_t>& map_ref_to_phys(
      SmallMat<space_dimT, n_vec, pt_coord_t>& points)
    {
      hy_assert(hyEdge_dimT == space_dimT, "This is only valid of the problem is of volumetype.");
      for (unsigned int i = 0; i < points.size(); ++i)
        hy_assert(points[i] >= 0. && points[i] <= 1., "Point must lie in reference square!");

      for (unsigned int j = 0; j < n_vec; ++j)
        for (unsigned int i = 0; i < space_dimT; ++i)
        {
          points(i, j) *= char_length[i];
          points(i, j) += translation[i];
        }

      return points;
    }
    /*!*********************************************************************************************
     * \brief   Return matrix column of the affine-linear transformation.
     *
     * \param   index   Index of the matrix column to be returned.
     * \retval  column  The specified matrix column.
     **********************************************************************************************/
    SmallVec<space_dimT, pt_coord_t> span_vec(const unsigned int index)
    {
      hy_assert(index < hyEdge_dimT, "There are only " << hyEdge_dimT << " spanning vectors.");

      SmallVec<space_dimT, pt_coord_t> span_vec;
      span_vec[dim_indices[index]] = char_length[index];

      return span_vec;
    }
    /*!*********************************************************************************************
     * \brief   Return reduced matrix R of the QR decomposition of the linear transoformation.
     **********************************************************************************************/
    const SmallSquareMat<hyEdge_dimT, pt_coord_t> mat_r()
    {
      return diagonal<hyEdge_dimT, hyEdge_dimT, pt_coord_t>(char_length);
    }
    /*!*********************************************************************************************
     * \brief   Return matrix Q of the QR decomposition of the linear transoformation.
     **********************************************************************************************/
    const SmallSquareMat<space_dimT, pt_coord_t> mat_q()
    {
      const auto& dim_ind = dim_indices.data();
      SmallSquareMat<space_dimT, pt_coord_t> mat_q;
      for (unsigned int index = 0; index < hyEdge_dimT; ++index)
        mat_q(dim_indices[index], index) = 1.;
      unsigned int helper = 0;
      for (unsigned int index = hyEdge_dimT; index < space_dimT; ++index)
      {
        while (std::find(dim_ind.begin(), dim_ind.end(), helper) != dim_ind.end())
          ++helper;
        mat_q(helper, index) = 1.;
        ++helper;
      }

      return mat_q;
    }
    /*!*********************************************************************************************
     * \brief   Return Haussdorff/Lebesque measure of the hyperedge.
     **********************************************************************************************/
    pt_coord_t area()
    {
      pt_coord_t area = 1.;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        area *= char_length[dim];
      return area;
    }
    /*!*********************************************************************************************
     * \brief   Return Haussdorff measure of the specified hypernode.
     **********************************************************************************************/
    pt_coord_t face_area(const unsigned int index)
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) faces.");
      pt_coord_t area = 1.;
      for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        if (dim != index / 2)
          area *= char_length[dim];
      return area;
    }
    /*!*********************************************************************************************
     * \brief   Return local normal of given index.
     *
     * Return outer unit normal with respect to the hypernode which is spanned by the vectors
     * spanning the phyiscal element, orthogonally projected to a hyEdge_dimT dimensional space,
     * but the vector of the given index. This is an element of the same dimension as the
     * reference element.
     **********************************************************************************************/
    Point<hyEdge_dimT, pt_coord_t> local_normal(const unsigned int index)
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) inner normals.");
      Point<hyEdge_dimT, pt_coord_t> normal;
      normal[index / 2] = index % 2 ? 1. : -1.;
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
    Point<space_dimT, pt_coord_t> inner_normal(const unsigned int index)
    {
      hy_assert(index < 2 * hyEdge_dimT, "A hyperedge has 2 * dim(hyEdge) inner normals.");
      Point<space_dimT, pt_coord_t> normal;
      normal[dim_indices[index / 2]] = index % 2 ? 1. : -1.;
      return normal;
    }
    /*!*********************************************************************************************
     * \brief   Return outer normal of given index.
     *
     * Return unit normal with respect to the hyperedge within the full space.
     **********************************************************************************************/
    Point<space_dimT, pt_coord_t> outer_normal(const unsigned int index)
    {
      hy_assert(index < space_dimT - hyEdge_dimT,
                "This function returns one of the dim(space) - dim(hyEdge) orthonormal vectors "
                  << "which are orthogonal to the hyperedge.");

      Point<space_dimT, pt_coord_t> outer_normal;
      unsigned int dim = 0;
      for (unsigned int cnt = 0; cnt < index + 1; ++dim)
        if (std::find(dim_indices.begin(), dim_indices.end(), dim) == dim_indices.end())
          ++cnt;
      outer_normal[dim] = 1.;

      return outer_normal;
    }
    /*!*********************************************************************************************
     * \brief   Return lexicographically ordered equidistant tensorial point of given index.
     **********************************************************************************************/
    template <unsigned int n_sub_points, typename one_dim_float_t>
    Point<space_dimT, pt_coord_t> lexicographic(
      unsigned int index,
      const SmallVec<n_sub_points, one_dim_float_t>& points_1d)
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
      return map_ref_to_phys(pt);
    }

    /*!*********************************************************************************************
     * \brief   Return equidistant tensorial point of given index on a given boundary face (slightly
     *          moved away from the boundary using boundary_scale), ordered lexicographically.
     **********************************************************************************************/
    template <unsigned int n_sub_points, typename one_dim_float_t>
    Point<space_dimT, pt_coord_t> boundary_lexicographic(
      unsigned int index,
      unsigned int boundary_number,
      float boundary_scale,
      const SmallVec<n_sub_points, one_dim_float_t>& points_1d)
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
      return map_ref_to_phys(pt);
    }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT   The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT       The dimension of the space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; }

 private:
  /*!***********************************************************************************************
   * \brief   Number of elements per spatial dimension.
   *
   * A vector / array comprising the number of elements in each spatial dimension.
   ************************************************************************************************/
  const SmallVec<space_dimT, unsigned int> n_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for elements.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>
    tpcc_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for refined elements.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>
    tpcc_ref_elem_;
  /*!***********************************************************************************************
   * \brief   Number of refined elements per element..
   ************************************************************************************************/
  const unsigned int n_loc_ref_elem;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class.
   *
   * The \c class \c UnitCube defines the geometry of the hypergraph. It contains the different
   * hyperedges (that actually are constructed everytime access is needed from e.g. the solver
   * class). Thus, its main purpose is to provide a structure that administrates the hyperedges that
   * are the return value of this structure.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
   *
   * To receive very general global loops, constructors need to account for the fact that the
   * specific topology / geometry of a hypergraph influences the way in which the hypergraph needs
   * to be constructed. The \c typedef implements the aspect, that a cubic hypergraph geometry which
   * is a unit cube is by default constructed by a vector / array that contains the numbers of
   * elements in the different dimensions.
   ************************************************************************************************/
  typedef ConstructorVecT constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a cubic that describes a unit cube hypergraph.
   *
   * Constructs a hypergraph from a \c constructor_value_type containing the elementens per spatial
   * dimension.
   *
   * \param   n_elements    The number of elements per spatial dimension.
   ************************************************************************************************/
  UnitCubeRefined(const constructor_value_type& n_elements)
  : n_elements_(n_elements),
    tpcc_elements_(
      Wrapper::create_tpcc<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        n_elements)),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervalsT))),
    n_loc_ref_elem(Hypercube<hyEdge_dimT>::pow(n_subintervalsT))
  {
  }
  /*!***********************************************************************************************
   * \brief   Construct a unit cube from its topological information
   *
   * Constructs a hypergraph from a \c Topology::Cubic containing the elementens per spatial
   * dimension.
   *
   * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
   ************************************************************************************************/
  UnitCubeRefined(const Topology::CubicRefined<hyEdge_dimT, space_dimT, n_subintervalsT>& other)
  : n_elements_(other.n_elements()),
    tpcc_elements_(
      Wrapper::create_tpcc<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        other.n_elements())),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervalsT))),
    n_loc_ref_elem(Hypercube<hyEdge_dimT>::pow(n_subintervalsT))
  {
  }
  /*!***********************************************************************************************
   * \brief   Get geometrical hyperedge of given index.
   *
   * This function returns the hyperedge of the given index, i.e., it returns the geometrical
   * hyperedge (\b not the topological information). The geometrical informatiom comprises the
   * indices of adjacent vertices (i.e. points) and information about their respective positions.
   *
   * \param   index       The index of the hyperedge to be returned.
   * \retval  hyperedge   Geometrical information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  const value_type operator[](const hyEdge_index_t index) const { return hyEdge(index, *this); }
};  // end class UnitCube

}  // end namespace Geometry
