#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/dense_la.hxx>

/*!*************************************************************************************************
 * \brief   A namespace containing classes describing hypergraph geometries.
 *
 * \todo    Discuss, whether we want to use array or SmallVec in internals.
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
template 
< 
  unsigned int hyEdge_dimT, unsigned int space_dimT, typename pt_coord_t = double,
  typename hyEdge_index_t = unsigned int
>
class UnitCube
{
  
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's hyperedge --- Edges of cubic HyperGraphs
   *          that form unit cubes.
   * 
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
    private:
      
      Point<space_dimT, pt_coord_t> translation;
      SmallVec<hyEdge_dimT, unsigned int> dim_indices;
      SmallVec<hyEdge_dimT, pt_coord_t> char_length;

      /*!*******************************************************************************************
       * \brief   Fill data of hyEdge.
       ********************************************************************************************/
      template<unsigned int hyEdge_dimTT>
      unsigned int fill_data
      ( 
        unsigned int index,
        const tpcc_elem_t<hyEdge_dimTT,space_dimT>& elem, const UnitCube& geometry
      )
      {
        if constexpr ( hyEdge_dimTT == 0 )
        {
          if (index == 0)
            for (unsigned int dim_pt = 0; dim_pt < space_dimT; ++dim_pt)
            {
              unsigned int ext_dim = exterior_direction<hyEdge_dimTT, space_dimT>(elem,dim_pt);
              translation[ext_dim] = (pt_coord_t) exterior_coordinate<hyEdge_dimTT, space_dimT>(elem, dim_pt)
                                        / (pt_coord_t) geometry.num_elements_[ext_dim];
              hy_assert( 0. <= translation[ext_dim] && translation[ext_dim] <= 1. ,
                          "The unit cube has only these cooridnates." );
            }
          unsigned int dim = 0;
          for (; dim < hyEdge_dimT && char_length[dim] != 0.; ++dim) ;
          if (index == (unsigned int) 1 << dim && char_length[dim] == 0.)
            for (unsigned int dim_pt = 0; dim_pt < space_dimT; ++dim_pt)
            {
              unsigned int ext_dim = exterior_direction<hyEdge_dimTT, space_dimT>(elem,dim_pt);
              pt_coord_t helper = (pt_coord_t) exterior_coordinate<hyEdge_dimTT, space_dimT>(elem, dim_pt)
                                        / (pt_coord_t) geometry.num_elements_[ext_dim];
              hy_assert( 0. <= helper && helper <= 1. ,
                          "The unit cube has only these cooridnates." );
              if ( helper != translation[ext_dim] )
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
          for (unsigned int i = 2 * hyEdge_dimTT - 2 ; i < 2 * hyEdge_dimTT; ++i)
            index = fill_data<hyEdge_dimTT-1>
                      ( index, get_face<hyEdge_dimTT, space_dimT>(elem, i), geometry );
        }
        return index;
      }
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
        fill_data<hyEdge_dimT>(0, elem, geometry);
      }
      
      /*!*******************************************************************************************
       * \brief   Map n_vec points from reference to physical element.
       ********************************************************************************************/
      template < unsigned int n_vec >
      SmallMat<space_dimT, n_vec, pt_coord_t> map_ref_to_phys
      (const SmallMat<hyEdge_dimT, n_vec, pt_coord_t>& pts)
      {
        for (unsigned int i = 0; i < pts.size(); ++i)
          hy_assert( pts[i] >= 0. && pts[i] <= 1. , "Point must lie in reference square!" );

        SmallMat<space_dimT, n_vec, pt_coord_t> phy_pts
          = rep_mat<space_dimT,n_vec,pt_coord_t>(translation);
        for (unsigned int j = 0; j < n_vec; ++j)
          for (unsigned int i = 0; i < hyEdge_dimT; ++i)
            phy_pts(dim_indices[i],j) += pts(i,j) * char_length[i];

        return phy_pts;
      }
      /*!*******************************************************************************************
       * \brief   Map n_vec points from reference to physical element.
       ********************************************************************************************/
      template < unsigned int n_vec >
      SmallMat<space_dimT, n_vec, pt_coord_t>& map_ref_to_phys
      (SmallMat<space_dimT, n_vec, pt_coord_t>& pts)
      {
        hy_assert( hyEdge_dimT == space_dimT ,
                   "This is only valid of the problem is of volumetype." )
        for (unsigned int i = 0; i < pts.size(); ++i)
          hy_assert( pts[i] >= 0. && pts[i] <= 1. , "Point must lie in reference square!");

        for (unsigned int j = 0; j < n_vec; ++j)
          for (unsigned int i = 0; i < space_dimT; ++i)
          {
            pts(i,j) *= char_length[i];
            pts(i,j) += translation[i];
          }

        return pts;
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

        SmallVec<space_dimT, pt_coord_t> span_vec;
        span_vec[dim_indices[index]] = char_length[index];

        return span_vec;
      }
      /*!*******************************************************************************************
       * \brief   Return reduced matrix R of the QR decomposition.
       ********************************************************************************************/
      const SmallSquareMat<hyEdge_dimT,pt_coord_t> mat_r()
      {
        return diagonal<hyEdge_dimT,hyEdge_dimT,pt_coord_t>(char_length);
      }
      /*!*******************************************************************************************
       * \brief   Return Haussdorff/Lebesque measure of the hyperedge.
       ********************************************************************************************/
      pt_coord_t area()
      {
        pt_coord_t area = 1.;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)  area *= char_length[dim];
        return area;
      }
      /*!*******************************************************************************************
       * \brief   Return Haussdorff measure of the specified hypernode.
       ********************************************************************************************/
      pt_coord_t face_area(const unsigned int index)
      {
        hy_assert( index < 2 * hyEdge_dimT ,
                   "A hyperedge has 2 * dim(hyEdge) faces." );
        pt_coord_t area = 1.;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          if (dim != index / 2)  area *= char_length[dim];
        return area;
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
        Point<hyEdge_dimT,pt_coord_t> normal;
        normal[index / 2] = index % 2 ? 1. : -1.;
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
        Point<space_dimT,pt_coord_t> normal;
        normal[dim_indices[index / 2]] = index % 2 ? 1. : -1.;
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

        Point<space_dimT,pt_coord_t> outer_normal;
        unsigned int dim = 0;
        for (unsigned int cnt = 0; cnt < index + 1; ++dim)
          if (std::find(dim_indices.begin(), dim_indices.end(), dim) == dim_indices.end()) ++cnt;
        outer_normal[dim] = 1.;

        return outer_normal;
      }
      /*!*******************************************************************************************
       * \brief   Return lexicographically ordered equidistant tensorial point of given index.
       ********************************************************************************************/
      template<unsigned int n_sub_points, typename one_dim_float_t>
      Point<space_dimT,pt_coord_t> lexicographic
      ( unsigned int index, const SmallVec<n_sub_points, one_dim_float_t>& points_1d )
      {
        static_assert( n_sub_points > 0 , "No subpoints do not make sense!" );
        hy_assert( index < std::pow(n_sub_points, hyEdge_dimT) ,
                   "The index must not exceed the number of prescribed lexicographic points." );
        Point<hyEdge_dimT,pt_coord_t> pt;
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          pt[dim] = (pt_coord_t) points_1d[index % n_sub_points];
          index /= n_sub_points;
        }
        return map_ref_to_phys(pt);
      }

    /*!*******************************************************************************************
     * \brief   Return equidistant tensorial point of given index on a given boundary (slightly moved
     * away from the boundary using boundary_scale), ordered lexicographically.
     ********************************************************************************************/
    template<unsigned int n_sub_points, typename one_dim_float_t>
    Point<space_dimT,pt_coord_t> boundary_lexicographic
        (unsigned int index, unsigned int boundary_number, float boundary_scale, const SmallVec<n_sub_points, one_dim_float_t>& points_1d )
    {
      static_assert( n_sub_points > 0 , "No subpoints do not make sense!" );
      hy_assert( index < std::pow(n_sub_points, hyEdge_dimT-1)*hyEdge_dimT*2 ,
                 "The index must not exceed the number of prescribed lexicographic points." );
      Point<hyEdge_dimT-1,pt_coord_t> subpt;
      index = index - std::pow(n_sub_points,hyEdge_dimT-1)*boundary_number;
      for (unsigned int subdim = 0; subdim < hyEdge_dimT-1; ++subdim)
      {
        subpt[subdim] = (pt_coord_t) points_1d[index % n_sub_points];
        index /= n_sub_points;
      }
      Point<hyEdge_dimT,pt_coord_t> pt;
      unsigned int subd = 0;
      for(unsigned int d = 0; d < hyEdge_dimT; d++)
      {
        if(boundary_number / 2 == d)
          pt[d] = boundary_scale*(boundary_number%2-0.5)+0.5;
        else
          pt[d] = subpt[subd++];
      }
      return map_ref_to_phys(pt);
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
    SmallVec<space_dimT, unsigned int> num_elements_;
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
    typedef SmallVec<space_dimT, unsigned int> constructor_value_type;
    /*!*********************************************************************************************
     * \brief   Construct a cubic that describes a cube hypergraph from a \c HyperGraph_Cubic.
     *
     * Constructs a hypergraph from a \c Topology::HyperGraph_Cubic containing the elementens per 
     * spatial dimension which is given as by its topology.
     * 
     * \param   other       The topology of the hypergraph that has the geometry of the unit cube.
     **********************************************************************************************/
    UnitCube(const constructor_value_type& num_elements)
    : num_elements_(num_elements), 
      tpcc_elements_(create_tpcc< hyEdge_dimT, space_dimT, hyEdge_index_t >(num_elements))
    { }
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
