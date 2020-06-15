# File containing geometric and topologic information for of Hypegraph being a triangle.

# The following parameters need to be given in the specified order!

Space_Dim     = 4;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).

N_Points      = 5;  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = 5;  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  = 10;  # Number of hyperedges that are contained in the hypergraph.


# After "POINTS:" the coordinates of the points have to appear. That is, the next lines conatain one
# point and therefore "Space_Dimension" numbers.

POINTS:
0.0 0.0 0.0 0.0
1.0 0.0 0.0 0.0
0.0 1.0 0.0 0.0
0.0 0.0 1.0 0.0
0.0 0.0 0.0 1.0

# After HYPERNODES_OF_HYPEREDGES the indices of hypernodes belonging to a hyperedge have to appear.

HYPERNODES_OF_HYPEREDGES:
0 1
0 2
0 3
0 4
1 2
1 3
1 4
2 3
2 4
3 4

# After TYPES OF HYPERFACES the type of the hyperface is denoted. As an example, 0 might indicate
# an interior face, 1 a Dirichlet face, and 2 a Neumann face.
# Note that a node can be an interior and a Neumann face, but if it is Dirichlet, all related faces
# have to be of Dirichlet type.

TYPES_OF_HYPERFACES:
0 1
0 2
0 3
0 4
1 2
1 3
1 4
2 3
2 4
3 4

# After POINTS_OF_HYPEREDGES the indices of points belonging to a hyperedge have to appear.

POINTS_OF_HYPEREDGES:
0 1
0 2
0 3
0 4
1 2
1 3
1 4
2 3
2 4
3 4