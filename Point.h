/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition of a class for point objects providing getter functions only.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


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
    Point();
    Point(const std::array<point_coord_type, space_dim>& coordinates);
    point_coord_type& operator[](const unsigned int coord_entry);
    bool operator==(const Point<space_dim>& other_point) const;
    bool operator!=(const Point<space_dim>& other_point) const;
    bool operator<(const Point<space_dim>& other_point) const;
};
/*
template<unsigned int space_dim>
double distance(const Point<space_dim>& left, const Point<space_dim>& right);
*/
#endif
