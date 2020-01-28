/*!*************************************************************************************************
 * \file    TypeDefs.hxx
 * \brief   Typedefs to abstract from standard C++ types.
 *
 * Defines several types, i.e. \c hypernode_index_type, \c hyperedge_index_type, \c dof_index_type,
 * and \c point_index_type which are supposed to be non-negative integers.
 * 
 * Additionally, \c dof_value_type and \c point_coord_type are floating types.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/

#ifndef TYPEDEFS_HXX
#define TYPEDEFS_HXX

typedef unsigned int  hypernode_index_type, hyperedge_index_type, dof_index_type, point_index_type;
typedef double        dof_value_type;
typedef float         point_coord_type;

#endif // end of ifndef TYPEDEFS_HXX
