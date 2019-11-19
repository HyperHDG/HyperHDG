/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "FuncAndQuad.h"
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;


template array<double, quadrature_points_amount(1)> quadrature_weights<1>();
template array<double, quadrature_points_amount(2)> quadrature_weights<2>();
template array<double, quadrature_points_amount(3)> quadrature_weights<3>();
template array<double, quadrature_points_amount(4)> quadrature_weights<4>();
template array<double, quadrature_points_amount(5)> quadrature_weights<5>();
template array<double, quadrature_points_amount(6)> quadrature_weights<6>();
template array<double, quadrature_points_amount(7)> quadrature_weights<7>();
template array<double, quadrature_points_amount(8)> quadrature_weights<8>();
template array<double, quadrature_points_amount(9)> quadrature_weights<9>();
template array<double, quadrature_points_amount(10)> quadrature_weights<10>();
template array<double, quadrature_points_amount(11)> quadrature_weights<11>();
template array<double, quadrature_points_amount(12)> quadrature_weights<12>();
template array<double, quadrature_points_amount(13)> quadrature_weights<13>();
template array<double, quadrature_points_amount(14)> quadrature_weights<14>();
template array<double, quadrature_points_amount(15)> quadrature_weights<15>();
template array<double, quadrature_points_amount(16)> quadrature_weights<16>();
template array<double, quadrature_points_amount(17)> quadrature_weights<17>();

template array< array<double, quadrature_points_amount(1)> , 1 + 1 > trial_functions_at_quadrature_points<1, 1>();
template array< array<double, quadrature_points_amount(2)> , 1 + 1 > trial_functions_at_quadrature_points<1, 2>();
template array< array<double, quadrature_points_amount(3)> , 1 + 1 > trial_functions_at_quadrature_points<1, 3>();
template array< array<double, quadrature_points_amount(4)> , 1 + 1 > trial_functions_at_quadrature_points<1, 4>();
template array< array<double, quadrature_points_amount(5)> , 1 + 1 > trial_functions_at_quadrature_points<1, 5>();
template array< array<double, quadrature_points_amount(6)> , 1 + 1 > trial_functions_at_quadrature_points<1, 6>();
template array< array<double, quadrature_points_amount(7)> , 1 + 1 > trial_functions_at_quadrature_points<1, 7>();
template array< array<double, quadrature_points_amount(8)> , 1 + 1 > trial_functions_at_quadrature_points<1, 8>();
template array< array<double, quadrature_points_amount(9)> , 1 + 1 > trial_functions_at_quadrature_points<1, 9>();
template array< array<double, quadrature_points_amount(10)> , 1 + 1 > trial_functions_at_quadrature_points<1, 10>();
template array< array<double, quadrature_points_amount(11)> , 1 + 1 > trial_functions_at_quadrature_points<1, 11>();
template array< array<double, quadrature_points_amount(12)> , 1 + 1 > trial_functions_at_quadrature_points<1, 12>();
template array< array<double, quadrature_points_amount(13)> , 1 + 1 > trial_functions_at_quadrature_points<1, 13>();
template array< array<double, quadrature_points_amount(14)> , 1 + 1 > trial_functions_at_quadrature_points<1, 14>();
template array< array<double, quadrature_points_amount(15)> , 1 + 1 > trial_functions_at_quadrature_points<1, 15>();
template array< array<double, quadrature_points_amount(16)> , 1 + 1 > trial_functions_at_quadrature_points<1, 16>();
template array< array<double, quadrature_points_amount(17)> , 1 + 1 > trial_functions_at_quadrature_points<1, 17>();
template array< array<double, quadrature_points_amount(1)> , 2 + 1 > trial_functions_at_quadrature_points<2, 1>();
template array< array<double, quadrature_points_amount(2)> , 2 + 1 > trial_functions_at_quadrature_points<2, 2>();
template array< array<double, quadrature_points_amount(3)> , 2 + 1 > trial_functions_at_quadrature_points<2, 3>();
template array< array<double, quadrature_points_amount(4)> , 2 + 1 > trial_functions_at_quadrature_points<2, 4>();
template array< array<double, quadrature_points_amount(5)> , 2 + 1 > trial_functions_at_quadrature_points<2, 5>();
template array< array<double, quadrature_points_amount(6)> , 2 + 1 > trial_functions_at_quadrature_points<2, 6>();
template array< array<double, quadrature_points_amount(7)> , 2 + 1 > trial_functions_at_quadrature_points<2, 7>();
template array< array<double, quadrature_points_amount(8)> , 2 + 1 > trial_functions_at_quadrature_points<2, 8>();
template array< array<double, quadrature_points_amount(9)> , 2 + 1 > trial_functions_at_quadrature_points<2, 9>();
template array< array<double, quadrature_points_amount(10)> , 2 + 1 > trial_functions_at_quadrature_points<2, 10>();
template array< array<double, quadrature_points_amount(11)> , 2 + 1 > trial_functions_at_quadrature_points<2, 11>();
template array< array<double, quadrature_points_amount(12)> , 2 + 1 > trial_functions_at_quadrature_points<2, 12>();
template array< array<double, quadrature_points_amount(13)> , 2 + 1 > trial_functions_at_quadrature_points<2, 13>();
template array< array<double, quadrature_points_amount(14)> , 2 + 1 > trial_functions_at_quadrature_points<2, 14>();
template array< array<double, quadrature_points_amount(15)> , 2 + 1 > trial_functions_at_quadrature_points<2, 15>();
template array< array<double, quadrature_points_amount(16)> , 2 + 1 > trial_functions_at_quadrature_points<2, 16>();
template array< array<double, quadrature_points_amount(17)> , 2 + 1 > trial_functions_at_quadrature_points<2, 17>();
template array< array<double, quadrature_points_amount(1)> , 3 + 1 > trial_functions_at_quadrature_points<3, 1>();
template array< array<double, quadrature_points_amount(2)> , 3 + 1 > trial_functions_at_quadrature_points<3, 2>();
template array< array<double, quadrature_points_amount(3)> , 3 + 1 > trial_functions_at_quadrature_points<3, 3>();
template array< array<double, quadrature_points_amount(4)> , 3 + 1 > trial_functions_at_quadrature_points<3, 4>();
template array< array<double, quadrature_points_amount(5)> , 3 + 1 > trial_functions_at_quadrature_points<3, 5>();
template array< array<double, quadrature_points_amount(6)> , 3 + 1 > trial_functions_at_quadrature_points<3, 6>();
template array< array<double, quadrature_points_amount(7)> , 3 + 1 > trial_functions_at_quadrature_points<3, 7>();
template array< array<double, quadrature_points_amount(8)> , 3 + 1 > trial_functions_at_quadrature_points<3, 8>();
template array< array<double, quadrature_points_amount(9)> , 3 + 1 > trial_functions_at_quadrature_points<3, 9>();
template array< array<double, quadrature_points_amount(10)> , 3 + 1 > trial_functions_at_quadrature_points<3, 10>();
template array< array<double, quadrature_points_amount(11)> , 3 + 1 > trial_functions_at_quadrature_points<3, 11>();
template array< array<double, quadrature_points_amount(12)> , 3 + 1 > trial_functions_at_quadrature_points<3, 12>();
template array< array<double, quadrature_points_amount(13)> , 3 + 1 > trial_functions_at_quadrature_points<3, 13>();
template array< array<double, quadrature_points_amount(14)> , 3 + 1 > trial_functions_at_quadrature_points<3, 14>();
template array< array<double, quadrature_points_amount(15)> , 3 + 1 > trial_functions_at_quadrature_points<3, 15>();
template array< array<double, quadrature_points_amount(16)> , 3 + 1 > trial_functions_at_quadrature_points<3, 16>();
template array< array<double, quadrature_points_amount(17)> , 3 + 1 > trial_functions_at_quadrature_points<3, 17>();
template array< array<double, quadrature_points_amount(1)> , 4 + 1 > trial_functions_at_quadrature_points<4, 1>();
template array< array<double, quadrature_points_amount(2)> , 4 + 1 > trial_functions_at_quadrature_points<4, 2>();
template array< array<double, quadrature_points_amount(3)> , 4 + 1 > trial_functions_at_quadrature_points<4, 3>();
template array< array<double, quadrature_points_amount(4)> , 4 + 1 > trial_functions_at_quadrature_points<4, 4>();
template array< array<double, quadrature_points_amount(5)> , 4 + 1 > trial_functions_at_quadrature_points<4, 5>();
template array< array<double, quadrature_points_amount(6)> , 4 + 1 > trial_functions_at_quadrature_points<4, 6>();
template array< array<double, quadrature_points_amount(7)> , 4 + 1 > trial_functions_at_quadrature_points<4, 7>();
template array< array<double, quadrature_points_amount(8)> , 4 + 1 > trial_functions_at_quadrature_points<4, 8>();
template array< array<double, quadrature_points_amount(9)> , 4 + 1 > trial_functions_at_quadrature_points<4, 9>();
template array< array<double, quadrature_points_amount(10)> , 4 + 1 > trial_functions_at_quadrature_points<4, 10>();
template array< array<double, quadrature_points_amount(11)> , 4 + 1 > trial_functions_at_quadrature_points<4, 11>();
template array< array<double, quadrature_points_amount(12)> , 4 + 1 > trial_functions_at_quadrature_points<4, 12>();
template array< array<double, quadrature_points_amount(13)> , 4 + 1 > trial_functions_at_quadrature_points<4, 13>();
template array< array<double, quadrature_points_amount(14)> , 4 + 1 > trial_functions_at_quadrature_points<4, 14>();
template array< array<double, quadrature_points_amount(15)> , 4 + 1 > trial_functions_at_quadrature_points<4, 15>();
template array< array<double, quadrature_points_amount(16)> , 4 + 1 > trial_functions_at_quadrature_points<4, 16>();
template array< array<double, quadrature_points_amount(17)> , 4 + 1 > trial_functions_at_quadrature_points<4, 17>();
template array< array<double, quadrature_points_amount(1)> , 5 + 1 > trial_functions_at_quadrature_points<5, 1>();
template array< array<double, quadrature_points_amount(2)> , 5 + 1 > trial_functions_at_quadrature_points<5, 2>();
template array< array<double, quadrature_points_amount(3)> , 5 + 1 > trial_functions_at_quadrature_points<5, 3>();
template array< array<double, quadrature_points_amount(4)> , 5 + 1 > trial_functions_at_quadrature_points<5, 4>();
template array< array<double, quadrature_points_amount(5)> , 5 + 1 > trial_functions_at_quadrature_points<5, 5>();
template array< array<double, quadrature_points_amount(6)> , 5 + 1 > trial_functions_at_quadrature_points<5, 6>();
template array< array<double, quadrature_points_amount(7)> , 5 + 1 > trial_functions_at_quadrature_points<5, 7>();
template array< array<double, quadrature_points_amount(8)> , 5 + 1 > trial_functions_at_quadrature_points<5, 8>();
template array< array<double, quadrature_points_amount(9)> , 5 + 1 > trial_functions_at_quadrature_points<5, 9>();
template array< array<double, quadrature_points_amount(10)> , 5 + 1 > trial_functions_at_quadrature_points<5, 10>();
template array< array<double, quadrature_points_amount(11)> , 5 + 1 > trial_functions_at_quadrature_points<5, 11>();
template array< array<double, quadrature_points_amount(12)> , 5 + 1 > trial_functions_at_quadrature_points<5, 12>();
template array< array<double, quadrature_points_amount(13)> , 5 + 1 > trial_functions_at_quadrature_points<5, 13>();
template array< array<double, quadrature_points_amount(14)> , 5 + 1 > trial_functions_at_quadrature_points<5, 14>();
template array< array<double, quadrature_points_amount(15)> , 5 + 1 > trial_functions_at_quadrature_points<5, 15>();
template array< array<double, quadrature_points_amount(16)> , 5 + 1 > trial_functions_at_quadrature_points<5, 16>();
template array< array<double, quadrature_points_amount(17)> , 5 + 1 > trial_functions_at_quadrature_points<5, 17>();

template array< array<double, quadrature_points_amount(1)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 1>();
template array< array<double, quadrature_points_amount(2)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 2>();
template array< array<double, quadrature_points_amount(3)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 3>();
template array< array<double, quadrature_points_amount(4)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 4>();
template array< array<double, quadrature_points_amount(5)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 5>();
template array< array<double, quadrature_points_amount(6)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 6>();
template array< array<double, quadrature_points_amount(7)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 7>();
template array< array<double, quadrature_points_amount(8)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 8>();
template array< array<double, quadrature_points_amount(9)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 9>();
template array< array<double, quadrature_points_amount(10)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 10>();
template array< array<double, quadrature_points_amount(11)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 11>();
template array< array<double, quadrature_points_amount(12)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 12>();
template array< array<double, quadrature_points_amount(13)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 13>();
template array< array<double, quadrature_points_amount(14)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 14>();
template array< array<double, quadrature_points_amount(15)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 15>();
template array< array<double, quadrature_points_amount(16)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 16>();
template array< array<double, quadrature_points_amount(17)> , 1 + 1 > derivs_of_trial_at_quadrature_points<1, 17>();
template array< array<double, quadrature_points_amount(1)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 1>();
template array< array<double, quadrature_points_amount(2)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 2>();
template array< array<double, quadrature_points_amount(3)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 3>();
template array< array<double, quadrature_points_amount(4)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 4>();
template array< array<double, quadrature_points_amount(5)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 5>();
template array< array<double, quadrature_points_amount(6)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 6>();
template array< array<double, quadrature_points_amount(7)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 7>();
template array< array<double, quadrature_points_amount(8)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 8>();
template array< array<double, quadrature_points_amount(9)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 9>();
template array< array<double, quadrature_points_amount(10)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 10>();
template array< array<double, quadrature_points_amount(11)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 11>();
template array< array<double, quadrature_points_amount(12)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 12>();
template array< array<double, quadrature_points_amount(13)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 13>();
template array< array<double, quadrature_points_amount(14)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 14>();
template array< array<double, quadrature_points_amount(15)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 15>();
template array< array<double, quadrature_points_amount(16)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 16>();
template array< array<double, quadrature_points_amount(17)> , 2 + 1 > derivs_of_trial_at_quadrature_points<2, 17>();
template array< array<double, quadrature_points_amount(1)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 1>();
template array< array<double, quadrature_points_amount(2)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 2>();
template array< array<double, quadrature_points_amount(3)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 3>();
template array< array<double, quadrature_points_amount(4)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 4>();
template array< array<double, quadrature_points_amount(5)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 5>();
template array< array<double, quadrature_points_amount(6)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 6>();
template array< array<double, quadrature_points_amount(7)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 7>();
template array< array<double, quadrature_points_amount(8)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 8>();
template array< array<double, quadrature_points_amount(9)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 9>();
template array< array<double, quadrature_points_amount(10)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 10>();
template array< array<double, quadrature_points_amount(11)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 11>();
template array< array<double, quadrature_points_amount(12)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 12>();
template array< array<double, quadrature_points_amount(13)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 13>();
template array< array<double, quadrature_points_amount(14)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 14>();
template array< array<double, quadrature_points_amount(15)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 15>();
template array< array<double, quadrature_points_amount(16)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 16>();
template array< array<double, quadrature_points_amount(17)> , 3 + 1 > derivs_of_trial_at_quadrature_points<3, 17>();
template array< array<double, quadrature_points_amount(1)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 1>();
template array< array<double, quadrature_points_amount(2)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 2>();
template array< array<double, quadrature_points_amount(3)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 3>();
template array< array<double, quadrature_points_amount(4)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 4>();
template array< array<double, quadrature_points_amount(5)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 5>();
template array< array<double, quadrature_points_amount(6)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 6>();
template array< array<double, quadrature_points_amount(7)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 7>();
template array< array<double, quadrature_points_amount(8)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 8>();
template array< array<double, quadrature_points_amount(9)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 9>();
template array< array<double, quadrature_points_amount(10)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 10>();
template array< array<double, quadrature_points_amount(11)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 11>();
template array< array<double, quadrature_points_amount(12)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 12>();
template array< array<double, quadrature_points_amount(13)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 13>();
template array< array<double, quadrature_points_amount(14)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 14>();
template array< array<double, quadrature_points_amount(15)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 15>();
template array< array<double, quadrature_points_amount(16)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 16>();
template array< array<double, quadrature_points_amount(17)> , 4 + 1 > derivs_of_trial_at_quadrature_points<4, 17>();
template array< array<double, quadrature_points_amount(1)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 1>();
template array< array<double, quadrature_points_amount(2)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 2>();
template array< array<double, quadrature_points_amount(3)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 3>();
template array< array<double, quadrature_points_amount(4)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 4>();
template array< array<double, quadrature_points_amount(5)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 5>();
template array< array<double, quadrature_points_amount(6)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 6>();
template array< array<double, quadrature_points_amount(7)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 7>();
template array< array<double, quadrature_points_amount(8)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 8>();
template array< array<double, quadrature_points_amount(9)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 9>();
template array< array<double, quadrature_points_amount(10)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 10>();
template array< array<double, quadrature_points_amount(11)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 11>();
template array< array<double, quadrature_points_amount(12)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 12>();
template array< array<double, quadrature_points_amount(13)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 13>();
template array< array<double, quadrature_points_amount(14)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 14>();
template array< array<double, quadrature_points_amount(15)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 15>();
template array< array<double, quadrature_points_amount(16)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 16>();
template array< array<double, quadrature_points_amount(17)> , 5 + 1 > derivs_of_trial_at_quadrature_points<5, 17>();

template array< array<double, 2> , 1 + 1 > trial_functions_at_boundaries<1>();
template array< array<double, 2> , 2 + 1 > trial_functions_at_boundaries<2>();
template array< array<double, 2> , 3 + 1 > trial_functions_at_boundaries<3>();
template array< array<double, 2> , 4 + 1 > trial_functions_at_boundaries<4>();
template array< array<double, 2> , 5 + 1 > trial_functions_at_boundaries<5>();

template array< array<double, 2> , 1 + 1 > derivs_of_trial_at_boundaries<1>();
template array< array<double, 2> , 2 + 1 > derivs_of_trial_at_boundaries<2>();
template array< array<double, 2> , 3 + 1 > derivs_of_trial_at_boundaries<3>();
template array< array<double, 2> , 4 + 1 > derivs_of_trial_at_boundaries<4>();
template array< array<double, 2> , 5 + 1 > derivs_of_trial_at_boundaries<5>();


double trial_function_eval(const unsigned int index, const double x_value)
{
  assert(0 <= index && index <= 5);
  assert(0. <= x_value && x_value <= 1.);
  switch (index)
  {
    case 0:  return 1.;
    case 1:  return sqrt(3) * (1. - 2. * x_value);
    case 2:  return sqrt(5) * ( (6. * x_value - 6.) * x_value + 1. );
    case 3:  return sqrt(7) * ( ( (20. * x_value - 30.) * x_value + 12. ) * x_value - 1. );
    case 4:  return sqrt(9) * ( ( ( (70. * x_value - 140.) * x_value + 90. ) * x_value - 20. ) * x_value + 1. );
    case 5:  return sqrt(11)* ( ( ( ( (252. * x_value - 630.) * x_value + 560. ) * x_value - 210. ) * x_value + 30. ) * x_value - 1. );
    default: assert(0 == 1);
  }
  assert(0 == 1);
  return 0.;
}


double deriv_of_trial_eval(const unsigned int index, const double x_value)
{
  assert(0 <= index && index <= 5);
  assert(0. <= x_value && x_value <= 1.);
  switch (index)
  {
    case 0:  return 0.;
    case 1:  return -sqrt(12);
    case 2:  return sqrt(5) * ( 12. * x_value - 6. );
    case 3:  return sqrt(7) * ( (60. * x_value - 60.) * x_value + 12. );
    case 4:  return sqrt(9) * ( ( (280. * x_value - 420.) * x_value + 180. ) * x_value - 20. );
    case 5:  return sqrt(11)* ( ( ( (1260. * x_value - 2520.) * x_value + 1680. ) * x_value - 420. ) * x_value + 30. );
    default: assert(0 == 1);
  }
  assert(0 == 1);
  return 0.;
}


template<unsigned int max_quad_degree>
array<double, quadrature_points_amount(max_quad_degree)> quadrature_points()
{
  constexpr unsigned int num_of_points = quadrature_points_amount(max_quad_degree);
  static_assert(1 <= num_of_points && num_of_points <= 9);
  array<double, num_of_points> quad_points;
  
  if constexpr (num_of_points == 1)
    quad_points = { 0. };
  if constexpr (num_of_points == 2)
    quad_points = { -sqrt(1./3.) , sqrt(1./3.) };
  if constexpr (num_of_points == 3)
    quad_points = { -sqrt(3./5.) , 0. , sqrt(3./5.) };
  if constexpr (num_of_points == 4)
    quad_points = { -sqrt(3./7.+2./7.*sqrt(6./5.)) , -sqrt(3./7.-2./7.*sqrt(6./5.)) ,
                     sqrt(3./7.-2./7.*sqrt(6./5.)) ,  sqrt(3./7.+2./7.*sqrt(6./5.)) };
  if constexpr (num_of_points == 5)
    quad_points = { -sqrt(5.+2.*sqrt(10./7.))/3. , -sqrt(5.-2.*sqrt(10./7.))/3., 0. ,
                     sqrt(5.-2.*sqrt(10./7.))/3. ,  sqrt(5.+2.*sqrt(10./7.))/3. };
  if constexpr (num_of_points == 6)
    quad_points = { 0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
                    0.2386191860831969, -0.9324695142031521,  0.9324695142031521 };
  if constexpr (num_of_points == 7)
    quad_points = { 0.0000000000000000,  0.4058451513773972, -0.4058451513773972,
                   -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,
                    0.9491079123427585 };
  if constexpr (num_of_points == 8)
    quad_points = {-0.1834346424956498,  0.1834346424956498, -0.5255324099163290,
                    0.5255324099163290, -0.7966664774136267,  0.7966664774136267,
                   -0.9602898564975363,  0.9602898564975363 };
  if constexpr (num_of_points == 9)
    quad_points = { 0.0000000000000000, -0.8360311073266358,  0.8360311073266358,
                   -0.9681602395076261,  0.9681602395076261, -0.3242534234038089,
                    0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
  /*                  
  switch (num_of_points)
  {
    case 1: quad_points = { 0. };
            break;
    case 2: quad_points = { -sqrt(1./3.) , sqrt(1./3.) };
            break;
    case 3: quad_points = { -sqrt(3./5.) , 0. , sqrt(3./5.) };
            break;
    case 4: quad_points = { -sqrt(3./7.+2./7.*sqrt(6./5.)) , -sqrt(3./7.-2./7.*sqrt(6./5.)) ,
                             sqrt(3./7.-2./7.*sqrt(6./5.)) ,  sqrt(3./7.+2./7.*sqrt(6./5.)) };
            break;
    case 5: quad_points = { -sqrt(5.+2.*sqrt(10./7.))/3. , -sqrt(5.-2.*sqrt(10./7.))/3., 0. ,
                             sqrt(5.-2.*sqrt(10./7.))/3. ,  sqrt(5.+2.*sqrt(10./7.))/3. };
            break;
    case 6: quad_points = { 0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
                            0.2386191860831969, -0.9324695142031521,  0.9324695142031521 };
            break;
    case 7: quad_points = { 0.0000000000000000,  0.4058451513773972, -0.4058451513773972,
                           -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,
                            0.9491079123427585 };
            break;
    case 8: quad_points = {-0.1834346424956498,  0.1834346424956498, -0.5255324099163290,
                            0.5255324099163290, -0.7966664774136267,  0.7966664774136267,
                           -0.9602898564975363,  0.9602898564975363 };
            break;
    case 9: quad_points = { 0.0000000000000000, -0.8360311073266358,  0.8360311073266358,
                           -0.9681602395076261,  0.9681602395076261, -0.3242534234038089,
                            0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
            break;
    default: assert( 0 == 1 );
  }
  */
  assert (num_of_points == quad_points.size() );
  // Transform quadrature points from [-1,1] -> [0,1]
  for_each(quad_points.begin(), quad_points.end(), [](double& pt){ pt = 0.5 * ( pt + 1. ); });
  return quad_points;
}

template<unsigned int max_quad_degree>
array<double, quadrature_points_amount(max_quad_degree)> quadrature_weights()
{
  constexpr unsigned int num_of_points = quadrature_points_amount(max_quad_degree);
  static_assert(1 <= num_of_points && num_of_points <= 9);
  array<double, num_of_points> quad_weights;
  
  if constexpr (num_of_points == 1)
    quad_weights = { 2. };
  if constexpr (num_of_points == 2)
    quad_weights = { 1. , 1. };
  if constexpr (num_of_points == 3)
    quad_weights = { 5./9. , 8./9. , 5./9. };
  if constexpr (num_of_points == 4)
    quad_weights = { 1./36.*(18. - sqrt(30.)) , 1./36.*(18. + sqrt(30.)) ,
                     1./36.*(18. + sqrt(30.)) , 1./36.*(18. - sqrt(30.)) };
  if constexpr (num_of_points == 5)
    quad_weights = { 1./900.*(322.-13.*sqrt(70.)) , 1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.+190.) ,
                     1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.-13.*sqrt(70.)) };
  if constexpr (num_of_points == 6)
    quad_weights = { 0.3607615730481386,  0.3607615730481386,  0.4679139345726910,
                     0.4679139345726910,  0.1713244923791704,  0.1713244923791700 };
  if constexpr (num_of_points == 7)
    quad_weights = { 0.4179591836734694,  0.3818300505051189,  0.3818300505051189,
                     0.2797053914892766,  0.2797053914892766,  0.1294849661688697,
                     0.1294849661688697 };
  if constexpr (num_of_points == 8)
    quad_weights = { 0.3626837833783620,  0.3626837833783620,  0.3137066458778873,
                     0.3137066458778873,  0.2223810344533745,  0.2223810344533745,
                     0.1012285362903763,  0.1012285362903763 };
  if constexpr (num_of_points == 9)
    quad_weights = { 0.3302393550012598,  0.1806481606948574,  0.1806481606948574,
                     0.0812743883615744,  0.0812743883615744,  0.3123470770400029,
                     0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
  /*  
  switch constexpr (num_of_points)
  {
    case 1: quad_weights = { 2. };
            break;
    case 2: quad_weights = { 1. , 1. };
            break;
    case 3: quad_weights = { 5./9. , 8./9. , 5./9. };
            break;
    case 4: quad_weights = { 1./36.*(18. - sqrt(30.)) , 1./36.*(18. + sqrt(30.)) ,
                             1./36.*(18. + sqrt(30.)) , 1./36.*(18. - sqrt(30.)) };
            break;
    case 5: quad_weights = { 1./900.*(322.-13.*sqrt(70.)) , 1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.+190.) ,
                             1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.-13.*sqrt(70.)) };
            break;
    case 6: quad_weights = { 0.3607615730481386,  0.3607615730481386,  0.4679139345726910,
                             0.4679139345726910,  0.1713244923791704,  0.1713244923791700 };
            break;
    case 7: quad_weights = { 0.4179591836734694,  0.3818300505051189,  0.3818300505051189,
                             0.2797053914892766,  0.2797053914892766,  0.1294849661688697,
                             0.1294849661688697 };
            break;
    case 8: quad_weights = { 0.3626837833783620,  0.3626837833783620,  0.3137066458778873,
                             0.3137066458778873,  0.2223810344533745,  0.2223810344533745,
                             0.1012285362903763,  0.1012285362903763 };
            break;
    case 9: quad_weights = { 0.3302393550012598,  0.1806481606948574,  0.1806481606948574,
                             0.0812743883615744,  0.0812743883615744,  0.3123470770400029,
                             0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
            break;
    default: assert( 0 == 1 );
  }
  */
  assert (num_of_points == quad_weights.size() );
  // Transform quadrature points from [-1,1] -> [0,1]
  for_each(quad_weights.begin(), quad_weights.end(), [](double& wt){ wt *= 0.5; });
  return quad_weights;
}


template<unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, quadrature_points_amount(max_quad_degree)> , max_poly_degree + 1 >
trial_functions_at_quadrature_points()
{
  constexpr unsigned int num_of_points = quadrature_points_amount(max_quad_degree);
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  static_assert( 1 <= num_of_points && num_of_points <= 9);
  
  array<double, num_of_points> quad_points = quadrature_points<max_quad_degree>();
  array< array<double, num_of_points> , max_poly_degree + 1 > fct_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
    for (unsigned int j = 0; j < num_of_points; ++j)
      fct_val[i][j] = trial_function_eval(i, quad_points[j]);
  
  return fct_val;
}


template<unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, quadrature_points_amount(max_quad_degree)> , max_poly_degree + 1 >
derivs_of_trial_at_quadrature_points()
{
  constexpr unsigned int num_of_points = quadrature_points_amount(max_quad_degree);
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  static_assert( 1 <= num_of_points && num_of_points <= 9);
  
  array<double, num_of_points> quad_points = quadrature_points<max_quad_degree>();
  array< array<double, num_of_points> , max_poly_degree + 1 > der_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
    for (unsigned int j = 0; j < num_of_points; ++j)
      der_val[i][j] = deriv_of_trial_eval(i, quad_points[j]);
  
  return der_val;
}


template<unsigned int max_poly_degree>
array< array<double, 2> , max_poly_degree + 1 > trial_functions_at_boundaries()
{
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  
  array< array<double, 2> , max_poly_degree + 1 > fct_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    fct_val[i][0] = trial_function_eval(i, 0.);
    fct_val[i][1] = trial_function_eval(i, 1.);
  }
  
  return fct_val;
}


template<unsigned int max_poly_degree>
array< array<double, 2> , max_poly_degree + 1 > derivs_of_trial_at_boundaries()
{
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  
  array< array<double, 2> , max_poly_degree + 1 > der_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    der_val[i][0] = deriv_of_trial_eval(i, 0.);
    der_val[i][1] = deriv_of_trial_eval(i, 1.);
  }
  
  return der_val;
}
