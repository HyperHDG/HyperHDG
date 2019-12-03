/*!*************************************************************************************************
 * @brief   This class implements a point in a d-dimensional space.
 * 
 * This class implements a point in a @f$d@f$-dimensional space, where the @f$d@f$ is given by the
 * template parameter @c space_dim.
 * 
 * @tparam  space_dim           The dimension of the space, the object is located in.
 * 
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/

#ifndef POINT_H
#define POINT_H

#include "TypeDefs.h"
#include <array>

template <unsigned int space_dim>
class Point
{
  private:
    std::array<point_coord_type, space_dim> coordinates_;
  public:
    /*!*********************************************************************************************
     * @brief   Empty constructor for a point.
     * 
     * Fills coordinates of the point with zeros.
     **********************************************************************************************/
    Point();
    /*!*********************************************************************************************
     * @brief   Construct point from array of coordinates.
     * 
     * Fills the point's array of coordinates with the input parameter. 
     * 
     * @param   coordinates   A @c std::array containing the coordinates of the point.
     **********************************************************************************************/
    Point(const std::array<point_coord_type, space_dim>& coordinates);
    /*!*********************************************************************************************
     * @brief   Return reference to single coordinate of a point.
     * 
     * @param   coord_entry   An @c unsigned @c int referring to the coordinate that is to be
     *                        returned.
     * @retval  coordinate    A reference to a @c point_coord_type describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    point_coord_type& operator[](const unsigned int coord_entry);
    /*!*********************************************************************************************
     * @brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns true if and only if both points
     * have exactly (that is not only with respect to some rounding errors) the same coordinates.
     * 
     * @param   other_point   Another @c Point<space_dim> that is to be dicriminate from.
     * @retval  isEqual       A @c boolean which is true if both points have the same coordinates.
     **********************************************************************************************/
    bool operator==(const Point<space_dim>& other_point) const;
    /*!*********************************************************************************************
     * @brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns false if and only if both
     * points have exactly (that is not only with respect to some rounding errors) the same
     * coordinates.
     * 
     * @param   other_point   Another @c Point<space_dim> that is to be dicriminate from.
     * @retval  isEqual       A @c boolean which is false if both points have the same coordinates.
     **********************************************************************************************/
    bool operator!=(const Point<space_dim>& other_point) const;
    /*!*********************************************************************************************
     * @brief   Find out whether the point is "smaller than" another point.
     * 
     * This function compares the point to another point and returns true if and only if the lowest
     * ranked coordinate (according to the coordinate index) where the both points are not equal of
     * the given point is smaller than that of the other point. It is false, if both points are
     * equal.
     * 
     * @param   other_point   Another @c Point<space_dim> that is to be dicriminate from.
     * @retval  isEqual       A @c boolean which is true if the left point is strictly smaller than
     *                        the right one.
     **********************************************************************************************/
    bool operator<(const Point<space_dim>& other_point) const;
};
/*
template<unsigned int space_dim>
double distance(const Point<space_dim>& left, const Point<space_dim>& right);
*/
#endif
