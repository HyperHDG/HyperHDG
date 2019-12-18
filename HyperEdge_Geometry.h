/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef HYPEREDGE_GEOMETRY_H
#define HYPEREDGE_GEOMETRY_H

#include "TypeDefs.h"
#include "Point.h"
#include <array>
//#include <vector>

namespace Geometry
{

template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperEdge_Cubic_UnitCube
{
  private:
//    std::array<hypernode_index_type, 2*hyperedge_dim> point_indices_;
    std::array<Point<space_dim>, 2*hyperedge_dim> points_;
//    std::array<bool, 2*hyperedge_dim> correct_hypernode_orientation_;
  public:
    HyperEdge_Cubic_UnitCube(const hyperedge_index_type index, const std::array<unsigned int, space_dim>& num_elements);
    Point<space_dim> point(unsigned int index) const;
//    std::vector<double> abs_det_of_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
//    std::vector< std::vector<double> > inv_of_transposed_jacobian_at_quad(const std::vector<double>& local_quadrature) const;
};

}

#endif
