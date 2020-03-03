/*!*************************************************************************************************
 * \file    examples_c++/PointTest.C
 * \brief   File that tests several aspects of the C++ implementation of Point (NOT FINISHED YET).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#include <HyperHDG/Point.hxx>

#include <array>
#include <random>

using namespace std;

/*!*************************************************************************************************
 * \brief   Function that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * \todo    Should we also add naive tests like checking whether return_zero_vector() returns vector
 *          of correct size only containing zeros?
 * 
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int space_dim, typename pt_coord_t> bool testPoint ( )
{
  static_assert( space_dim != 0 , "Space dimension must be strictly larger than zero!" );
  bool success = true;
  
  const pt_coord_t minR = -10.;
  const pt_coord_t maxR = +10.;

  Point<space_dim,pt_coord_t> ptA, ptB, pt;
  pt_coord_t randNrC, result, helper;
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<pt_coord_t> dis(minR, maxR);
  
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == 0. && ptB[dim] == 0. ,
               "Initially, all points are supposed to be filled with zeros only!" );
    if (ptA[dim] != 0. || ptB[dim] != 0.)  success = false;
    do    ptA[dim] = dis(gen);
    while ( ptA[dim] == 0. );
    do    ptB[dim] = dis(gen);
    while ( ptB[dim] == 0. );
  }
  
  do  randNrC = dis(gen);
  while ( randNrC == 0. );

  const pt_coord_t randNrCC = randNrC;
  const Point<space_dim,pt_coord_t> ptAC(ptA), ptBC(ptB);

  hy_assert( ptA == ptAC && !(ptA != ptAC) && !(ptA < ptAC) && !(ptAC < ptA) ,
             "Points created from same data should be equal even if one is const and one is not! "
             << "This implies that they are neither unequal nor smaller than one another." );
  if (ptA != ptAC || !(ptA == ptAC) || ptA < ptAC || ptAC < ptA)  success = false;

  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] ,
               "Operator[] must return same values for Point and const Point if they are equal!" );
    if (ptA[dim] != ptAC[dim])  success = false;
  }
  
  Point<space_dim,pt_coord_t> ptMC( move(ptA) );
  Point<space_dim,pt_coord_t> ptMA = Point<space_dim,pt_coord_t>(ptAC);
  ptA = ptAC;
  
  hy_assert( ptMC == ptMA && ptMA == ptA ,
             "Move constructor and move assignment should give the same points!" );
  if (ptMC != ptMA || ptMA != ptA)  success = false;
  
  if constexpr (space_dim == 1)       ptA[0] -= 1e-5;
  else if constexpr (space_dim == 2)  ptA[1] -= 1e-5;
  else                              { ptA[1] -= 1e-5; ptA[2] += 1.; }

  hy_assert( !(ptA == ptAC) && ptA != ptAC && ptA < ptAC && !(ptAC < ptA) ,
             "Manipulated point is smaller than old point and therefore unequal!" );
  if (!(ptA != ptAC) || ptA == ptAC || !(ptA < ptAC) || ptAC < ptA)  success = false;

  hy_assert( ptA != ptAC && ptA != ptMC && ptA != ptMA && ptAC == ptMC && ptAC == ptMA ,
             "Other constructors should have created independent points!" )
  if (ptA == ptAC || ptA == ptMC || ptA == ptMA || ptAC != ptMC || ptAC != ptMA)  success = false;

  ptA = ptAC;
  
  
  ptA += randNrCC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] + randNrCC , "Operator+= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] + randNrCC)  success = false;
  }
  ptA = ptAC;
  
  ptA -= randNrCC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] - randNrCC , "Operator-= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] - randNrCC)  success = false;
  }
  ptA = ptAC;
  
  ptA *= randNrCC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] * randNrCC , "Operator*= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] * randNrCC)  success = false;
  }
  ptA = ptAC;
  
  ptA /= randNrCC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] / randNrCC , "Operator/= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] / randNrCC)  success = false;
  }
  ptA = ptAC;
  
  
  ptA += ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] + ptBC[dim] , "Operator+= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] + ptBC[dim])  success = false;
  }
  ptA = ptAC;
  
  ptA -= ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] - ptBC[dim] , "Operator-= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] - ptBC[dim])  success = false;
  }
  ptA = ptAC;
  
  ptA *= ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] * ptBC[dim] , "Operator*= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] * ptBC[dim])  success = false;
  }
  ptA = ptAC;
  
  ptA /= ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptA[dim] == ptAC[dim] / ptBC[dim] , "Operator/= failed specified behaviour!" );
    if (ptA[dim] != ptAC[dim] / ptBC[dim])  success = false;
  }
  ptA = ptAC;
  

  result = ptAC * ptBC;
  helper = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)  helper += ptA[dim] * ptB[dim];
  hy_assert( result == helper , "Operator* (scalar product) failed specified behaviour!" );
  if (result != helper)  success = false;
  
  
  pt = ptAC + ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] + ptBC[dim] , "Operator+ failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] + ptBC[dim])  success = false;
  }
  
  pt = ptAC - ptBC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] - ptBC[dim] , "Operator- failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] - ptBC[dim])  success = false;
  }
  
  pt = hada_prod(ptAC, ptBC);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] * ptBC[dim] , "hada_prod failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] * ptBC[dim])  success = false;
  }
  
  pt = hada_divi(ptAC, ptBC);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] / ptBC[dim] , "hada_divi failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] / ptBC[dim])  success = false;
  }
  

  pt = randNrC + ptAC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC + ptAC[dim] , "Operator+ failed specified behaviour!" );
    if (pt[dim] != randNrC + ptAC[dim])  success = false;
  }
  
  pt = ptAC + randNrC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC + ptAC[dim] , "Operator+ failed specified behaviour!" );
    if (pt[dim] != randNrC + ptAC[dim])  success = false;
  }
  
  pt = randNrC - ptAC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC - ptAC[dim] , "Operator- failed specified behaviour!" );
    if (pt[dim] != randNrC - ptAC[dim])  success = false;
  }

  pt = ptAC - randNrC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] - randNrC , "Operator- failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] - randNrC)  success = false;
  }

  pt = randNrC * ptAC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC * ptAC[dim] , "Operator* failed specified behaviour!" );
    if (pt[dim] != randNrC * ptAC[dim])  success = false;
  }

  pt = ptAC * randNrC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC * ptAC[dim] , "Operator* failed specified behaviour!" );
    if (pt[dim] != randNrC * ptAC[dim])  success = false;
  }
   
  pt = randNrC / ptAC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == randNrC / ptAC[dim] , "Operator/ failed specified behaviour!" );
    if (pt[dim] != randNrC / ptAC[dim])  success = false;
  }

  pt = ptAC / randNrC;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( pt[dim] == ptAC[dim] / randNrC , "Operator/ failed specified behaviour!" );
    if (pt[dim] != ptAC[dim] / randNrC)  success = false;
  }
  
  
  hy_assert( norm_1(ptAC) == norm_p(ptAC, (pt_coord_t) 1.) 
               && norm_2(ptAC) == norm_p(ptAC, (pt_coord_t) 2.) ,
             "Norms should be the same indpendent of their implementation!" );
  if(norm_1(ptAC) != norm_p(ptAC, (pt_coord_t) 1.) || norm_2(ptAC) != norm_p(ptAC, (pt_coord_t) 2.))
    success = false;

  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    hy_assert( ptAC[dim] <= norm_infty(ptAC) && ptAC[dim] >= -norm_infty(ptAC) ,
               "Each component (in absolute value) should be smaller than infty norm." );
    if (ptAC[dim] > norm_infty(ptAC) || ptAC[dim] < -norm_infty(ptAC))  success = false;
  }
//  hy_assert( norm_infty(ptAC) < norm_p(ptAC, 10) ,
//             "Infty norm is smallest of all norms!" );
//  if (norm_infty(ptAC) > norm_p(ptAC, 1000))  success = false;
    
  return success;
}


int main(int argc, char *argv[])
{
  bool success = true, worked;
  
  worked = testPoint<1,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 1 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<2,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 2 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<3,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 3 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<4,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 4 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<5,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 5 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<6,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 6 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<7,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 7 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<8,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 8 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<9,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 9 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<10,double>();
  hy_assert( worked , "Testing point implemantation with double precision of dimension 10 failed!");
  if (!worked)  success = false;
  
  
  worked = testPoint<1,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 1 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<2,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 2 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<3,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 3 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<4,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 4 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<5,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 5 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<6,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 6 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<7,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 7 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<8,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 8 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<9,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 9 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<10,float>();
  hy_assert( worked , "Testing point implemantation with float precision of dimension 10 failed!" );
  if (!worked)  success = false;

  return success - 1;
}
