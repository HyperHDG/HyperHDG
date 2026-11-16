#pragma once  // Ensure that file is included only once in a single compilation.

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Kind tags: compile-time selection of which condensed operator a local solver builds.
 *
 * A local solver provides one overload of \c local_operator(edge, Kind, args) per Kind it supports
 * and, if it supports others than \c Stiffness, lists them as \c Kinds (a std::tuple of the tags;
 * a kind is any default-constructible type). The global loop takes the kind as a template argument
 * of prepare and the sweeps, Stiffness by default, and forwards the tag unchanged. Only
 * \c Stiffness exists in this prototype.
 **************************************************************************************************/
struct Stiffness
{
};
}  // end of namespace GlobalLoop
