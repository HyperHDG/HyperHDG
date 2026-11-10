#pragma once  // Ensure that file is included only once in a single compilation.

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Kind tags: compile-time selection of which condensed operator a local solver builds.
 *
 * A local solver provides one overload of \c local_operator(edge, Kind, args) per Kind it supports.
 * The global loop forwards the tag unchanged. Only \c Stiffness exists in this prototype.
 **************************************************************************************************/
struct Stiffness
{
};
}  // end of namespace GlobalLoop
