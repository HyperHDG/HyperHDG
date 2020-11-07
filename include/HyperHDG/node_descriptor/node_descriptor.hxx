#pragma once  // Ensure that file is included only once in a single compilation.

// This file contains doxygen information only!

/*!*************************************************************************************************
 * \brief   A namespace containing different classes describing the types of hypernodes.
 *
 * A hypergraph consists of hyperedges (where the PDE is usually defined) and hypernodes. Usually,
 * these nodes encode the continuity of primal and dual variables. However, nodes might also of
 * some boundary type or have additional properties.
 *
 * To discriminate between several types of hypernodes, the node descriptors return the type of each
 * node belonging to a hyperedge (similarly to the topology which returns the hypernode indices).
 * The local solver knowns, which type (from the node descriptor) encodes which property (such as a
 * Dirichlet node, node where certain continuity conditions are supposed to be met, etc.).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
namespace NodeDescriptor
{
}  // end of namespace NodeDescriptor
