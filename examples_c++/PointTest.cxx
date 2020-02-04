/*!*************************************************************************************************
 * \file    examples_c++/PointTest.C
 * \brief   File that tests several aspects of the C++ implementation of Point (NOT FINISHED YET).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#include <Point.hxx>

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
template<unsigned int space_dim> bool testPoint ( )
{
  static_assert( space_dim != 0 , "Space dimension must be strictly larger than zero!" );
  bool successful = true;
  
  const unsigned int array_len = 5;
  const pt_coord_t minR = -10.;
  const pt_coord_t maxR = +10.;
  
  pt_coord_t randNr, result, aux;
  Point<space_dim> helper, res;
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<pt_coord_t> dis(minR, maxR);
  
  array< array<pt_coord_t, space_dim>, array_len > coord_array;
  array< Point<space_dim>, array_len > pt_array;
  
  for (unsigned int i = 0; i < array_len; ++i)
  {
    for (unsigned int j = 0; j < space_dim; ++j)
    {
      do  coord_array[i][j] = dis(gen);
      while ( coord_array[i][j] == 0. );
      hy_assert( pt_array[i][j] == 0. ,
                 "Initially, all points are supposed to be filled with zeros only!" );
      if (pt_array[i][j] != 0.)  successful = false;
    }
    pt_array[i] = coord_array[i];
    for (unsigned int dim = 0; dim < space_dim; ++dim)
      hy_assert( pt_array[0][dim] == coord_array[0][dim] && pt_array[0][dim] != 0. ,
                 "Point created from array should contain array's coordinates not equal to zero." );
  }
  
  do  randNr = dis(gen);
  while ( randNr == 0. );
  const pt_coord_t randNrC = randNr;
               
  Point<space_dim> pt(pt_array[0]);
  hy_assert( pt == pt_array[0] && !(pt != pt_array[0]) && !(pt < pt_array[0]) ,
             "Point should be equal to point it is created from and not unequal to that point!" );
  
  const Point<space_dim> ptC(coord_array[0]);
  hy_assert( pt == ptC ,
             "Points created from same data should be equal even if one is const and one is not!" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] ,
               "Operator[] should return same values for Point and const Point!" );
  
  Point<space_dim> ptMC( move(pt) );
  Point<space_dim> ptMA = Point<space_dim>(coord_array[0]);
  pt = ptC;
  
  hy_assert( ptMC == ptMA && ptMA == pt ,
             "Move constructor and move assignment should give the same points!" );
  
  if constexpr (space_dim == 1)       pt[0] -= 1e-5;
  else if constexpr (space_dim == 2)  pt[1] -= 1e-5;
  else                              { pt[1] -= 1e-5; pt[2] += 1.; }
  hy_assert( !(pt == pt_array[0]) && pt != pt_array[0] && pt < pt_array[0] ,
             "Manipulated point is smaller than old point and therefor unequal!" );
  hy_assert( pt != ptC && pt != ptMC && pt != ptMA && ptC == ptMC && ptC == ptMA ,
             "Other constructors should have created independent points!" )
  pt = ptC;
  
  
  pt += randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] + randNr ,
               "Operator+= should add same number to all components!" );
  pt = ptC;
  
  pt -= randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] - randNr ,
               "Operator-= should subtract same number from all components!" );
  pt = ptC;
  
  pt *= randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] * randNr ,
               "Operator*= should multiply same number to all components!" );
  pt = ptC;
  
  pt /= randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] / randNr ,
               "Operator/= divide all components by randNr!" );
  pt = ptC;
  
  
  helper = pt_array[1];
  pt += helper;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] + pt_array[1][dim] ,
               "Operator+= should add up two points without changing the argument!" );
  hy_assert( helper == pt_array[1] ,
             "Operator+= should add up two points without changing the argument!" );
  pt = ptC;
  
  helper = pt_array[2];
  pt -= helper;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] - pt_array[2][dim] ,
               "Operator-= should subtract argument without changing it!" );
  hy_assert( helper == pt_array[2] ,
             "Operator+= should add up two points without changing the argument!" );
  pt = ptC;
  
  helper = pt_array[3];
  pt *= helper;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] * pt_array[3][dim] ,
               "Operator*= should multiply with the argument without changing it!" );
  hy_assert( helper == pt_array[3] ,
             "Operator+= should add up two points without changing the argument!" );
  pt = ptC;
  
  helper = pt_array[4];
  pt /= helper;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == ptC[dim] / pt_array[4][dim] ,
               "Operator/= divide by the argument without changing it!" );
  hy_assert( helper == pt_array[4] ,
             "Operator+= should add up two points without changing the argument!" );
  pt = ptC;
  
  helper = pt_array[1];
  result = pt * pt_array[1];
  aux = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)  aux += pt[dim] * helper[dim];
  hy_assert( result == aux && helper == pt_array[1] && pt == ptC ,
             "Operator* should determine the scalar product without changing the arguments!" );
  
  
  pt = pt_array[0] + pt_array[1];
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == pt_array[0][dim] + pt_array[1][dim] ,
               "At least one of the operator+ failed!" <<  pt[dim] << "  " << pt_array[0][dim] + (pt_array[1][dim] + pt_array[2][dim]) + pt_array[3][dim] + (pt_array[4][dim] + pt_array[5][dim]));
  
  pt = pt_array[0] - pt_array[1];
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == pt_array[0][dim] - pt_array[1][dim] ,
               "At least one of the operator- failed!" );
  
  pt = hada_prod(pt_array[0], pt_array[1]);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == pt_array[0][dim] * pt_array[1][dim] ,
               "At least one of the hada_prod failed!" );
  
  pt = hada_divi(pt_array[0], pt_array[1]);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( pt[dim] == pt_array[0][dim] / pt_array[1][dim] ,
               "At least one of the hada_divi failed!" );
  
  pt = ptC;
  res = randNr + pt;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr + pt[dim] , "One of the scalar multiplications failed!" );
  
  pt = ptC;
  res = pt + randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr + pt[dim] , "One of the scalar multiplications failed!" );
   
  
  pt = ptC;
  res = randNr - pt;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr - pt[dim] , "One of the scalar multiplications failed!" );
  
  
  pt = ptC;
  res = pt - randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == pt[dim] - randNr , "One of the scalar multiplications failed!" );
  
  
  pt = ptC;
  res = randNr * pt;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr * pt[dim] , "One of the scalar multiplications failed!" );
  
  pt = ptC;
  res = pt * randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr * pt[dim] , "One of the scalar multiplications failed!" );
   
  
  pt = ptC;
  res = randNr / pt;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == randNr / pt[dim] , "One of the scalar multiplications failed!" );
  
  
  pt = ptC;
  res = pt / randNr;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    hy_assert( res[dim] == pt[dim] / randNr , "One of the scalar multiplications failed!" );
  
  // TODO: Norms and stream
  
    
  return successful;
}


int main(int argc, char *argv[])
{
  bool success = true, worked;
  
  worked = testPoint<1>();
  hy_assert( worked , "Testing point implemantation of dimension 1 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<2>();
  hy_assert( worked , "Testing point implemantation of dimension 2 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<3>();
  hy_assert( worked , "Testing point implemantation of dimension 3 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<4>();
  hy_assert( worked , "Testing point implemantation of dimension 4 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<5>();
  hy_assert( worked , "Testing point implemantation of dimension 5 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<6>();
  hy_assert( worked , "Testing point implemantation of dimension 6 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<7>();
  hy_assert( worked , "Testing point implemantation of dimension 7 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<8>();
  hy_assert( worked , "Testing point implemantation of dimension 8 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<9>();
  hy_assert( worked , "Testing point implemantation of dimension 9 failed!" );
  if (!worked)  success = false;
  
  worked = testPoint<10>();
  hy_assert( worked , "Testing point implemantation of dimension 10 failed!" );
  if (!worked)  success = false;
  
  
  return !success;
}
