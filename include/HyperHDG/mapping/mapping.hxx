#pragma once  // Ensure that file is included only once in a single compilation.

// This file contains doxygen information only!

/*!*************************************************************************************************
 * \brief   Namespace for mappings from reference elements to physical elements, etc.
 *
 * Mappings of a \c hyEdge_dimT dimensional unit square to a subset of a subset of the real numbers
 * to the power \c space_dimT heavily reliy on QR decomposition. Additionally, there are four
 * different types of normals:
 *
 * - The normals of the unit square in \c hyEdge_dimT dimensions.
 * - The normals of R times the unit square in \c hyEdge_dimT dimensions (called local_normal).
 * - The normals of QR times the unit square within the planar spanned by the columns of the
 *   transformation matrix (denoted inner normals) in \c space_dimT dimensions.
 * - The orthonormal vectors to the planar spanned by the transformation matrix (outer normals).
 *
 * The QR decomposition has a special normalization, i.e., the matrix Q suffices det(Q) = +1, i.e.
 * Q describes a movement (no mirrioring), and R has non-negative diagonal entries --- except for
 * the entry (0,0) which may have negative sign. Thus, the sign of entry (0,0) describes, whether
 * matrix a is orientation preserving or not.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
namespace Mapping
{
}  // end of namespace Mapping
