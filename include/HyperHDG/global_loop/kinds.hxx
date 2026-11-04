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
/*!*************************************************************************************************
 * \brief   Tag for an operator that is used once (a sweep without prepare) instead of stored.
 *
 * A local solver may provide local_operator(edge, Kind, args, Streaming) for a representation that
 * is cheap to build (e.g. LU factors instead of the condensed matrix); without it, the loop uses
 * local_operator(edge, Kind, args) for streaming as well.
 **************************************************************************************************/
struct Streaming
{
};
}  // end of namespace GlobalLoop
