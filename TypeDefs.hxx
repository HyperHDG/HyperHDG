/*!*************************************************************************************************
 * \file    TypeDefs.hxx
 * \brief   Typedefs to abstract from standard C++ types.
 *
 * Defines several types, i.e. \c hyNode_index_t, \c hyEdge_index_t, \c dof_index_type, and 
 * \c pt_index_t which are supposed to be non-negative integers.
 * 
 * Additionally, \c dof_value_t and \c pt_coord_t are floating types.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/

#ifndef TYPEDEFS_HXX
#define TYPEDEFS_HXX

typedef unsigned int  hyNode_index_t, hyEdge_index_t, dof_index_type, pt_index_t;
typedef double        dof_value_t;

#endif // end of ifndef TYPEDEFS_HXX
