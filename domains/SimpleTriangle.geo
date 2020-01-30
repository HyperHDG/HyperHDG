# File containing geometric and topologic information for of Hypegraph being a triangle.

# The following parameters need to be given in the specified order!

Space_Dim     = 2;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).

N_Points      = 3;  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = 3;  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  = 3;  # Number of hyperedges that are contained in the hypergraph.


# After "POINTS:" the coordinates of the points have to appear. That is, the next lines conatain one
# point and therefore "Space_Dimension" numbers.

POINTS:
0.0 0.0
2.0 0.0
1.0 1.0

# After "HYPERNODES_OF_HYPEREDGES the indices of hypernodes belonging to a hyperedge have to appear.

HYPERNODES_OF_HYPEREDGES:
0 1
0 2
1 2

# After "POINTS_OF_HYPEREDGES: the indices of points belonging to a hyperedge have to appear.

POINTS_OF_HYPEREDGES:
0 1
0 2
1 2
