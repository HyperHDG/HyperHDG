#pragma once  // Ensure that file is included only once in a single compilation.

// This file contains doxygen information only!

/*!*************************************************************************************************
 * \brief   A namespace containing classes describing hypergraph geometries.
 *
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain \f$\Omega\f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * There are four different types of normals specified (cf. the respective mapping):
 * - The normals of the reference element.
 * - The normals of of the transformed element which is still an element of the hyEdge_dimT
 *   dimensional space, but already has the shape of the physical element (called local_normal).
 * - The space_dimT dimensional normals of the phyiscal element located within the planar (or curve
 *   if curved elements are allowed) which is spanned by the element.
 * - The orthonormal vectors to the planar/curve spanned by phyiscal element (outer normals).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
namespace Geometry
{
}  // end namespace Geometry
