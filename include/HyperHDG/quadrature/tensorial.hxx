#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/quadrature/one_dimensional.hxx>

#include <array>
#include <cmath>
#include <numeric>

namespace Quadrature
{
/*!*************************************************************************************************
 * \brief   General integrator class on tensorial hypergraphs.
 *
 * \tparam  max_poly_degree   Maximum degree of evaluated polynomials.
 * \tparam  max_quad_degree   Desired degree of accuracy.
 * \tparam  quadrature_t      The quadrature rule applied.
 * \param   shape_t           Type of one-dimensional shape functions.
 * \tparam  return_t          Floating type specification. Default is double.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename quadrature_t, typename shape_t, typename return_t = double>
struct Tensorial
{
  typedef shape_t shape_fun_t;
  static constexpr unsigned int dim() { return shape_t::dim(); }
  static constexpr unsigned int n_fun_1D = shape_t::shape_fun_t::shape_fun_1d::n_fun();
  /*!***********************************************************************************************
   * \brief   Calculate the amount of quadrature points.
   *
   * \tparam  quadrature_t      The quadrature rule applied.
   * \param   max_quad_degree   Desired degree of accuracy.
   * \param   local_dimensions  Dimension of the underlying domain. Defaullt is one.
   * \retval  n_quad_points     Amount of needed quadrature points.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  static constexpr unsigned int n_quad_points()
  {
    return Hypercube<dim()>::pow(quadrature_t::n_quad_points());
  }

  /*!***********************************************************************************************
   * \brief   Quadrature points on one-dimensional unit interval.
   *
   * Returns the quadrature points of the quadrature rule with accuracy order \c max_quad_degree on
   * a one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  max_quad_degree   Desired degree of accuracy.
   * \tparam  quadrature_t      The quadrature rule applied.
   * \tparam  return_t          Floating type specification. Default is double.
   * \retval  quad_points       \c std::array containing the quadrature points.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  static inline std::array<return_t, quadrature_t::n_points()> quad_points()
  {
    return quadrature_t::template points<return_t>();
  }
  /*!***********************************************************************************************
   * \brief   Quadrature weights on one-dimensional unit interval.
   *
   * Returns the quadrature weights of the quadrature rule with accuracy order \c max_quad_degree on
   * a one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  max_quad_degree   Desired degree of accuracy.
   * \tparam  quadrature_t      The quadrature rule applied.
   * \tparam  return_t          Floating type specification. Default is double.
   * \retval  quad_weights      \c std::array containing the quadrature weights.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  static inline std::array<return_t, quadrature_t::n_points()> quad_weights()
  {
    return quadrature_t::template weights<return_t>();
  }

  // Shape functions & their derivatives evaluated at quadrature's points:

  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at quadrature points.
   *
   * Returns the values of the shape functions on \f$[0,1]\f$ of degree at most
   * \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
   * one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  max_poly_degree   Maximum degree of evaluated polynomials.
   * \tparam  max_quad_degree   Desired degree of accuracy.
   * \tparam  quadrature_t      The quadrature rule applied.
   * \tparam   shape_t           Type of one-dimensional shape functions.
   * \tparam  return_t          Floating type specification. Default is double.
   * \retval  quad_vals         \c std::array of polynomial degrees containing \c std::array of
   *                            quadrature points (the shape functions are evaluated at).
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  static inline std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D>
  shape_fcts_at_quad_points()
  {
    std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < quadrature_t::n_points(); ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template fct_val<return_t>(
          fct, quadrature_t::points()[pt]);

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at quadrature points.
   *
   * Returns the values of the derivatives of orthonormal shape functions on \f$[0,1]\f$ of degree
   *at most \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
   * one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  max_poly_degree   Maximum degree of evaluated polynomials.
   * \tparam  max_quad_degree   Desired degree of accuracy.
   * \tparam  quadrature_t      The quadrature rule applied.
   * \tparam   shape_t           Type of one-dimensional shape functions.
   * \tparam  return_t          Floating type specification. Default is double.
   * \retval  quad_vals         \c std::array of polynomial degrees containing \c std::array of
   *                            quadrature points (the shape functions' derivatives are evaluated).
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  std::array<std::array<return_t, quadrature_t::n_points()>,
             n_fun_1D> static inline shape_ders_at_quad_points()
  {
    std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < quadrature_t::n_points(); ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template der_val<return_t>(
          fct, quadrature_t::points()[pt]);

    return result;
  }

  std::array<std::array<return_t, 2>, n_fun_1D> static inline shape_fcts_at_bdrs()
  {
    std::array<std::array<return_t, 2>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < 2; ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template fct_val<return_t>(fct, pt);

    return result;
  }

  std::array<std::array<return_t, 2>, n_fun_1D> static inline shape_ders_at_bdrs()
  {
    std::array<std::array<return_t, 2>, n_fun_1D> result;

    for (unsigned int fct = 0; fct < n_fun_1D; ++fct)
      for (unsigned int pt = 0; pt < 2; ++pt)
        result[fct][pt] = shape_t::shape_fun_t::shape_fun_1d::template der_val<return_t>(fct, pt);

    return result;
  }

  const std::array<return_t, quadrature_t::n_points()> quad_points_, quad_weights_;
  const std::array<std::array<return_t, quadrature_t::n_points()>, n_fun_1D> shape_fcts_at_quad_,
    shape_ders_at_quad_;
  const std::array<std::array<return_t, 2>, n_fun_1D> trial_bdr_, trial_der_bdr_;
  /*!***********************************************************************************************
   * \brief   Constructor for general integrator class.
   ************************************************************************************************/
  Tensorial()
  : quad_points_(quad_points()),
    quad_weights_(quad_weights()),
    shape_fcts_at_quad_(shape_fcts_at_quad_points()),
    shape_ders_at_quad_(shape_ders_at_quad_points()),
    trial_bdr_(shape_fcts_at_bdrs()),
    trial_der_bdr_(shape_ders_at_bdrs())
  {
    hy_assert(quad_weights_.size() == quad_points_.size(),
              "Number of quadrature weights and points must be equal!");
    hy_assert(shape_fcts_at_quad_.size() == shape_ders_at_quad_.size(),
              "Number of shape functions and their derivatives must be equal!");
    for (unsigned int i = 0; i < shape_fcts_at_quad_.size(); ++i)
      hy_assert(quad_points_.size() == shape_fcts_at_quad_[i].size() &&
                  shape_fcts_at_quad_[i].size() == shape_ders_at_quad_[i].size(),
                "Number of quadrature points needs to be equal in all cases!");
  }

  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape functions.
   *
   * \param   i             Local index of local one-dimensional shape function.
   * \param   j             Local index of local one-dimensional shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  inline return_t integrate_1D_phiphi(const unsigned int i, const unsigned int j) const
  {
    hy_assert(i < shape_fcts_at_quad_.size() && j < shape_fcts_at_quad_.size(),
              "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights_.size(); ++q)
      result += quad_weights_[q] * shape_fcts_at_quad_[i][q] * shape_fcts_at_quad_[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape function and one derivative.
   *
   * \param   i             Local index of local one-dimensional shape function.
   * \param   j             Local index of local one-dimensional shape function (with derivative).
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  inline return_t integrate_1D_phiDphi(const unsigned int i, const unsigned int j) const
  {
    hy_assert(i < shape_fcts_at_quad_.size() && j < shape_fcts_at_quad_.size(),
              "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights_.size(); ++q)
      result += quad_weights_[q] * shape_fcts_at_quad_[i][q] * shape_ders_at_quad_[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of one-dimensional shape function and one derivative.
   *
   * \param   i             Local index of local one-dimensional shape function (with derivative).
   * \param   j             Local index of local one-dimensional shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  inline return_t integrate_1D_Dphiphi(const unsigned int i, const unsigned int j) const
  {
    hy_assert(i < shape_fcts_at_quad_.size() && j < shape_fcts_at_quad_.size(),
              "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights_.size(); ++q)
      result += quad_weights_[q] * shape_ders_at_quad_[i][q] * shape_fcts_at_quad_[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of two one-dimensional shape functions' derivatives.
   *
   * \param   i             Local index of local one-dimensional shape function (with derivative).
   * \param   j             Local index of local one-dimensional shape function (with derivative).
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  inline return_t integrate_1D_DphiDphi(const unsigned int i, const unsigned int j) const
  {
    hy_assert(i < shape_fcts_at_quad_.size() && j < shape_fcts_at_quad_.size(),
              "Indices of shape functions must be smaller than amount of shape functions.");
    return_t result = 0.;

    for (unsigned int q = 0; q < quad_weights_.size(); ++q)
      result += quad_weights_[q] * shape_ders_at_quad_[i][q] * shape_ders_at_quad_[j][q];

    return result;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  return_t integrate_vol_phiphi(const unsigned int i, const unsigned int j) const
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function amd derivative over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function (with derivative).
   * \param   dim_der       Dimension of the derivative.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  return_t integrate_vol_phiDphi(const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int dim_der) const
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (dim_der == dim_fct)
        integral *= integrate_1D_phiDphi(dec_i[dim_fct], dec_j[dim_fct]);
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function and derivative over dimT-dimensional unit volume.
   *
   * \param   i             Local index of local shape function (with derivative).
   * \param   j             Local index of local shape function.
   * \param   dim_der       Dimension of the derivative.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  return_t integrate_vol_Dphiphi(const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int dim_der) const
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (dim_der == dim_fct)
        integral *= integrate_1D_Dphiphi(dec_i[dim_fct], dec_j[dim_fct]);
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over dimT-dimensional volume's boundary.
   *
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  return_t integrate_bdr_phiphi(const unsigned int i,
                                const unsigned int j,
                                const unsigned int bdr) const
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, dim()> dec_j = Hypercube<dim()>::index_decompose(j, n_fun_1D);
    unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (bdr_dim == dim_fct)
        integral *= trial_bdr_[dec_i[dim_fct]][bdr_ind] * trial_bdr_[dec_j[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function of volume times shape function of volume's face
   *          over dimT-dimensional volume's boundary.
   *
   * \param   i             Local index of local volume shape function.
   * \param   j             Local index of local boundary shape function.
   * \param   bdr           Boundary face index.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  return_t integrate_bdr_phipsi(const unsigned int i,
                                const unsigned int j,
                                const unsigned int bdr) const
  {
    return_t integral = 1.;
    std::array<unsigned int, dim()> dec_i = Hypercube<dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(dim() - 1, 1U)> dec_j =
      Hypercube<dim() - 1>::index_decompose(j, n_fun_1D);
    unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < dim(); ++dim_fct)
      if (bdr_dim == dim_fct)
        integral *= trial_bdr_[dec_i[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct - (dim_fct > bdr_dim)]);
    return integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  func          Function that is also to be integrated.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   time          Time at which function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_vol_phiphifunc(const unsigned int i,
                                    const unsigned int j,
                                    GeomT& geom,
                                    const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_val *= quad_weights_[dec_q[dim]] * shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]] *
                    shape_fcts_at_quad_[dec_j[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  func          Function that is also to be integrated.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   dimension     Local dimension with respect to which the vector-function is integrated.
   * \param   geom          Geometrical information.
   * \param   time          Time at which function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            SmallVec<GeomT::space_dim(), return_t> fun(const Point<GeomT::space_dim(), return_t>&,
                                                       const return_t)>
  return_t integrate_vol_phiphivecfunc(const unsigned int i,
                                       const unsigned int j,
                                       const unsigned int dimension,
                                       GeomT& geom,
                                       const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;
    const SmallSquareMat<GeomT::space_dim(), return_t> mat_q =
      (SmallSquareMat<GeomT::space_dim(), return_t>)geom.mat_q();

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_val *= quad_weights_[dec_q[dim]] * shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]] *
                    shape_fcts_at_quad_[dec_j[dim]][dec_q[dim]];
      }
      integral +=
        scalar_product(fun(geom.map_ref_to_phys(quad_pt), time), mat_q.get_column(dimension)) *
        quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT>
  return_t integrate_vol_phiphi(const unsigned int i, const unsigned int j, GeomT& geom) const
  {
    constexpr unsigned int dimT = GeomT::hyEdge_dim();
    return_t integral = 1.;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    for (unsigned int dim_fct = 0; dim_fct < dimT; ++dim_fct)
      integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of linear combinations of shape functions over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  array_size    Size of arrays containing coefficients of linear combinations.
   * \tparam  floating_t    The floating point type for the calculation.
   * \param   is            Coefficients of local shape functions.
   * \param   js            Coefficients of local shape functions.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of lineat combinations of shape functions.
   ************************************************************************************************/
  template <typename GeomT, std::size_t array_size, typename floating_t>
  return_t integrate_vol_phiphi(const std::array<floating_t, array_size>& is,
                                const std::array<floating_t, array_size>& js,
                                GeomT& geom) const
  {
    return_t integral = 0., quad_val, is_val, js_val, val_helper;

    std::array<unsigned int, GeomT::hyEdge_dim()> dec_k, dec_q;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      is_val = 0.;
      js_val = 0.;

      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
        quad_val *= quad_weights_[dec_q[dim]];

      for (unsigned int k = 0; k < array_size; ++k)
      {
        dec_k = Hypercube<GeomT::hyEdge_dim()>::index_decompose(k, n_fun_1D);
        val_helper = 1.;
        for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
          val_helper *= shape_fcts_at_quad_[dec_k[dim]][dec_q[dim]];
        is_val += is[k] * val_helper;
        js_val += js[k] * val_helper;
      }
      integral += quad_val * is_val * js_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape function times some function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Function whose product with shape function is integrated.
   * \param   i             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   time          Time at which function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_vol_phifunc(const unsigned int i, GeomT& geom, const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_val *= quad_weights_[dec_q[dim]] * shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral * geom.area();
  }
  /*!***********************************************************************************************
   * \brief   Average integral of product of shape function times some function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Function whose product with shape function is integrated.
   * \param   i             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   time          Time at which function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_volUni_phifunc(const unsigned int i,
                                    GeomT& geom,
                                    const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, quadrature_t::n_points());
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_val *= quad_weights_[dec_q[dim]] * shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]];
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times other shape function over some geometry.
   *
   * \note    poly_deg_i and poly_deg_j must be set to max_poly_degree (which is also their default)
   *          if the basis is not hierarchic.
   * \todo    Generalize and recover correct implementation!
   *
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  poly_deg_i    Polynomial degree of shape functions associated to i.
   * \tparam  poly_deg_j    Polynomial degree of shape functions associated to j.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            unsigned int poly_deg_i = shape_t::degree(),
            unsigned int poly_deg_j = shape_t::degree()>
  SmallVec<GeomT::hyEdge_dim(), return_t> integrate_vol_nablaphiphi(const unsigned int i,
                                                                    const unsigned int j,
                                                                    GeomT& geom) const
  {
    static_assert(poly_deg_i <= shape_t::degree() && poly_deg_j <= shape_t::degree(),
                  "The maximum polynomial degrees must be larger than or equal to the given ones.");
    SmallVec<GeomT::hyEdge_dim(), return_t> integral(1.);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, poly_deg_i + 1);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, poly_deg_j + 1);
    for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
        if (dim == dim_fct)
          integral[dim] *= integrate_1D_Dphiphi(dec_i[dim_fct], dec_j[dim_fct]);
        else
          integral[dim] *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return geom.area() * integral / transposed(geom.mat_r());
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape functions times other function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   geom          Geometrical information.
   * \param   time          Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_vol_nablaphinablaphifunc(const unsigned int i,
                                              const unsigned int j,
                                              GeomT& geom,
                                              return_t time = 0.) const
  {
    return_t integral = 0., quad_weight;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt, nabla_phi_i, nabla_phi_j;
    const SmallMat<GeomT::hyEdge_dim(), GeomT::hyEdge_dim(), return_t> rrT =
      mat_times_transposed_mat(geom.mat_r(), geom.mat_r());

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      nabla_phi_i = 1.;
      nabla_phi_j = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_weight *= quad_weights_[dec_q[dim]];
        for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct)
          {
            nabla_phi_i[dim_fct] *= shape_ders_at_quad_[dec_i[dim]][dec_q[dim]];
            nabla_phi_j[dim_fct] *= shape_ders_at_quad_[dec_j[dim]][dec_q[dim]];
          }
          else
          {
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]];
            nabla_phi_j[dim_fct] *= shape_fcts_at_quad_[dec_j[dim]][dec_q[dim]];
          }
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), time) *
                  scalar_product(nabla_phi_i, nabla_phi_j / rrT);
    }
    return geom.area() * integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate derivative of shape function times other function over some geometry.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \param   i             Local index of local shape function.
   * \param   dim_der       Dimension with respect to which derivative is calculated.
   * \param   geom          Geometrical information.
   * \param   time          Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_vol_derphifunc(const unsigned int i,
                                    const unsigned int dim_der,
                                    GeomT& geom,
                                    return_t time = 0.) const
  {
    return_t integral = 0., quad_weight;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt, nabla_phi_i;
    const SmallSquareMat<GeomT::hyEdge_dim(), return_t> rT = transposed(geom.mat_r());

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_weight *= quad_weights_[dec_q[dim]];
        for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad_[dec_i[dim]][dec_q[dim]];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]];
        }
      }
      integral +=
        quad_weight * fun(geom.map_ref_to_phys(quad_pt), time) * (nabla_phi_i / rT)[dim_der];
    }
    return geom.area() * integral;
  }
  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times shape function times other function times
   *          normal over some geometry's boundary.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   bdr           Index of the boundatry face to integrate over.
   * \param   geom          Geometrical information.
   * \param   time          Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_bdr_nablaphiphinufunc(const unsigned int i,
                                           const unsigned int j,
                                           const unsigned int bdr,
                                           GeomT& geom,
                                           return_t time = 0.) const
  {
    return_t integral = 0., quad_weight, phi_j;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, std::max(GeomT::hyEdge_dim() - 1, 1U)> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt, nabla_phi_i, normal;
    const SmallMat<GeomT::hyEdge_dim(), GeomT::hyEdge_dim(), return_t> rT =
      transposed(geom.mat_r());
    const unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      phi_j = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        if (dim == bdr_dim)
        {
          normal[dim] = 2. * bdr_ind - 1.;
          quad_pt[dim] = (return_t)bdr_ind;
          phi_j *= trial_bdr_[dec_j[dim]][bdr_ind];
        }
        else
        {
          normal[dim] = 0.;
          quad_pt[dim] = quad_points_[dec_q[dim - (dim > bdr_dim)]];
          phi_j *= shape_fcts_at_quad_[dec_j[dim]][dec_q[dim - (dim > bdr_dim)]];
          quad_weight *= quad_weights_[dec_q[dim - (dim > bdr_dim)]];
        }
        for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct && dim == bdr_dim)
            nabla_phi_i[dim_fct] *= trial_der_bdr_[dec_i[dim]][bdr_ind];
          else if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad_[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
          else if (dim == bdr_dim)
            nabla_phi_i[dim_fct] *= trial_bdr_[dec_i[dim]][bdr_ind];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad_[dec_i[dim]][dec_q[dim - (dim > bdr)]];
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), time) * phi_j *
                  scalar_product(normal, nabla_phi_i / rT);
    }
    return geom.face_area(bdr) * integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate gradient of shape function times shape function times other function times
   *          normal over some geometry's boundary.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  fun           Weight function that is additionally integrated.
   * \param   i             Local index of local shape function with gradient.
   * \param   j             Local index of local shape function.
   * \param   bdr           Index of the boundatry face to integrate over.
   * \param   geom          Geometrical information.
   * \param   time          Time at which fun is evaluated.
   * \retval  integral      Integral of product of both shape function's weighted gradients.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_bdr_nablaphipsinufunc(const unsigned int i,
                                           const unsigned int j,
                                           const unsigned int bdr,
                                           GeomT& geom,
                                           return_t time = 0.) const
  {
    return_t integral = 0., quad_weight, phi_j;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(GeomT::hyEdge_dim() - 1, 1U)> dec_j =
      Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(j, n_fun_1D);
    std::array<unsigned int, std::max(GeomT::hyEdge_dim() - 1, 1U)> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt, nabla_phi_i, normal;
    const SmallMat<GeomT::hyEdge_dim(), GeomT::hyEdge_dim(), return_t> rT =
      transposed(geom.mat_r());
    const unsigned int bdr_dim = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      phi_j = 1.;
      nabla_phi_i = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        if (dim == bdr_dim)
        {
          normal[dim] = 2. * bdr_ind - 1.;
          quad_pt[dim] = (return_t)bdr_ind;
        }
        else
        {
          normal[dim] = 0.;
          quad_pt[dim] = quad_points_[dec_q[dim - (dim > bdr_dim)]];
          phi_j *= shape_fcts_at_quad_[dec_j[dim - (dim > bdr_dim)]][dec_q[dim - (dim > bdr_dim)]];
          quad_weight *= quad_weights_[dec_q[dim - (dim > bdr_dim)]];
        }
        for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
        {
          if (dim == dim_fct && dim == bdr_dim)
            nabla_phi_i[dim_fct] *= trial_der_bdr_[dec_i[dim]][bdr_ind];
          else if (dim == dim_fct)
            nabla_phi_i[dim_fct] *= shape_ders_at_quad_[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
          else if (dim == bdr_dim)
            nabla_phi_i[dim_fct] *= trial_bdr_[dec_i[dim]][bdr_ind];
          else
            nabla_phi_i[dim_fct] *= shape_fcts_at_quad_[dec_i[dim]][dec_q[dim - (dim > bdr_dim)]];
        }
      }
      integral += quad_weight * fun(geom.map_ref_to_phys(quad_pt), time) * phi_j *
                  scalar_product(normal, nabla_phi_i / rT);
    }
    return geom.face_area(bdr) * integral;
  }

  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions over boundary face.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \param   i             Local index of local shape function.
   * \param   j             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT>
  return_t integrate_bdr_phiphi(const unsigned int i,
                                const unsigned int j,
                                const unsigned int bdr,
                                GeomT& geom) const
  {
    return_t integral = 1.;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_j =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(j, n_fun_1D);
    unsigned int dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
      if (dim == dim_fct)
        integral *= trial_bdr_[dec_i[dim_fct]][bdr_ind] * trial_bdr_[dec_j[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct]);
    return integral * geom.face_area(bdr);
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions of volumen and skeletal over boundary face.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \param   i             Local index of local volumne shape function.
   * \param   j             Local index of local skeletal shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT>
  return_t integrate_bdr_phipsi(const unsigned int i,
                                const unsigned int j,
                                const unsigned int bdr,
                                GeomT& geom) const
  {
    return_t integral = 1.;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(GeomT::hyEdge_dim() - 1, 1U)> dec_j =
      Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(j, n_fun_1D);
    unsigned int dim = bdr / 2, bdr_ind = bdr % 2;
    for (unsigned int dim_fct = 0; dim_fct < GeomT::hyEdge_dim(); ++dim_fct)
      if (dim == dim_fct)
        integral *= trial_bdr_[dec_i[dim_fct]][bdr_ind];
      else
        integral *= integrate_1D_phiphi(dec_i[dim_fct], dec_j[dim_fct - (dim_fct > dim)]);
    return integral * geom.face_area(bdr);
  }
  /*!***********************************************************************************************
   * \brief   Integrate product of shape functions times some function over boundary face.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  func          Function that is multiplied by shape function.
   * \param   i             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \param   time          Time at which the function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_bdr_phifunc(const unsigned int i,
                                 const unsigned int bdr,
                                 GeomT& geom,
                                 const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i =
      Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
    std::array<unsigned int, std::max(1U, GeomT::hyEdge_dim() - 1)> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;
    unsigned int dim_bdr = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(q, n_fun_1D);
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        if (dim == dim_bdr)
        {
          quad_pt[dim] = bdr_ind;
          quad_val *= trial_bdr_[dec_i[dim]][bdr_ind];
        }
        else
        {
          quad_pt[dim] = quad_points_[dec_q[dim - (dim > dim_bdr)]];
          quad_val *= quad_weights_[dec_q[dim - (dim > dim_bdr)]] *
                      shape_fcts_at_quad_[dec_i[dim]][dec_q[dim - (dim > dim_bdr)]];
        }
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral * geom.face_area(bdr);
  }

  /*template
  < typename GeomT, return_t fun(const Point<GeomT::space_dim(),return_t>&, const return_t) >
  return_t integrate_bdrUni_phifunc
  (const unsigned int i, const unsigned int bdr, GeomT& geom, const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i = index_decompose<GeomT::hyEdge_dim()>(i);
    std::array<unsigned int, std::max(1U,GeomT::hyEdge_dim()-1)> dec_q;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;
    unsigned int dim_bdr = bdr / 2 , bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()-1); ++q)
    {
      dec_q = index_decompose
                <GeomT::hyEdge_dim()-1,quadrature_t::compute_n_quad_points(max_quad_degree)>(q);
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        if (dim == dim_bdr)
        {
          quad_pt[dim] = bdr_ind;
          quad_val *= trial_bdr_[dec_i[dim]][bdr_ind];
        }
        else
        {
          quad_pt[dim] = quad_points_[dec_q[dim - (dim > dim_bdr)]];
          quad_val *= quad_weights_[dec_q[dim - (dim > dim_bdr)]]
                        * shape_fcts_at_quad_[dec_i[dim]][dec_q[dim - (dim > dim_bdr)]];
        }
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral;
  }*/
  /*!***********************************************************************************************
   * \brief   Average integral of product of skeletal shape functions times some function.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  func          Function that is multiplied by shape function.
   * \param   i             Local index of local shape function.
   * \param   bdr           Boundary face index.
   * \param   geom          Geometrical information.
   * \param   time          Time at which the function is evaluated.
   * \retval  integral      Integral of product of both shape functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t)>
  return_t integrate_bdrUni_psifunc(const unsigned int i,
                                    const unsigned int bdr,
                                    GeomT& geom,
                                    const return_t time = 0.) const
  {
    return_t integral = 0., quad_val;
    std::array<unsigned int, std::max(1U, GeomT::hyEdge_dim() - 1)> dec_q,
      dec_i = Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(i, n_fun_1D);
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;
    unsigned int dim_bdr = bdr / 2, bdr_ind = bdr % 2;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim() - 1); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim() - 1>::index_decompose(q, n_fun_1D);
      quad_val = 1.;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        if (dim == dim_bdr)
          quad_pt[dim] = bdr_ind;
        else
        {
          quad_pt[dim] = quad_points_[dec_q[dim - (dim > dim_bdr)]];
          quad_val *=
            quad_weights_[dec_q[dim - (dim > dim_bdr)]] *
            shape_fcts_at_quad_[dec_i[dim - (dim > dim_bdr)]][dec_q[dim - (dim > dim_bdr)]];
        }
      }
      integral += fun(geom.map_ref_to_phys(quad_pt), time) * quad_val;
    }
    return integral;
  }
  /*!***********************************************************************************************
   * \brief   Squared L2 distance of some function and an discrete function on volume.
   *
   * \tparam  GeomT         Geometry which is the integration domain.
   * \tparam  func          Function whose distance is measured.
   * \param   coeffs        Coefficients of discrete function.
   * \param   geom          Geometrical information.
   * \param   time          Time at which the function is evaluated.
   * \retval  integral      Squared distance of functions.
   ************************************************************************************************/
  template <typename GeomT,
            return_t fun(const Point<GeomT::space_dim(), return_t>&, const return_t),
            std::size_t n_coeff>
  return_t integrate_vol_diffsquare_discana(const std::array<return_t, n_coeff> coeffs,
                                            GeomT& geom,
                                            const return_t time = 0.) const
  {
    return_t integral = 0., quad_weight;
    std::array<unsigned int, GeomT::hyEdge_dim()> dec_i, dec_q;
    std::array<return_t, n_coeff> quad_val;
    Point<GeomT::hyEdge_dim(), return_t> quad_pt;

    for (unsigned int q = 0; q < std::pow(quad_weights_.size(), GeomT::hyEdge_dim()); ++q)
    {
      dec_q = Hypercube<GeomT::hyEdge_dim()>::index_decompose(q, n_fun_1D);
      quad_weight = 1.;
      quad_val = coeffs;
      for (unsigned int dim = 0; dim < GeomT::hyEdge_dim(); ++dim)
      {
        quad_pt[dim] = quad_points_[dec_q[dim]];
        quad_weight *= quad_weights_[dec_q[dim]];
        for (unsigned int i = 0; i < n_coeff; ++i)
        {
          dec_i = Hypercube<GeomT::hyEdge_dim()>::index_decompose(i, n_fun_1D);
          quad_val[i] *= shape_fcts_at_quad_[dec_i[dim]][dec_q[dim]];
        }
      }
      integral += quad_weight * std::pow(fun(geom.map_ref_to_phys(quad_pt), time) -
                                           std::accumulate(quad_val.begin(), quad_val.end(), 0.),
                                         2);
    }
    return integral * geom.area();
  }
};  // end of class Integrator

}  // end of namespace Quadrature