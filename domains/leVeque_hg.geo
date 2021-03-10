# File containing geometric and topologic information for of Hypegraph being a triangle.

# The following parameters need to be given in the specified order!

Space_Dim     = 3;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 2;  # Dimension of hyperedge (must be uniform).

N_Points      = 12;  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = 18;  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  =  6;  # Number of hyperedges that are contained in the hypergraph.


# After "POINTS:" the coordinates of the points have to appear. That is, the next lines conatain one
# point and therefore "Space_Dimension" numbers.

POINTS:
-2.0   0.0   0.0
-1.0   0.0   0.0
-1.0   0.0   1.0
-2.0   0.0   1.0
 0.0   1.0   0.0
 0.0   1.0   1.0
 0.0  -1.0   0.0
 0.0  -1.0   1.0
 1.0   0.0   0.0
 2.0   0.0   0.0
 2.0   0.0   1.0
 1.0   0.0   1.0


# After HYPERNODES_OF_HYPEREDGES the indices of hypernodes belonging to a hyperedge have to appear.

HYPERNODES_OF_HYPEREDGES:
 0  1  2  3
 1  4  5  6
 1  7  8  9
 4 10 11 12
 7 10 13 14
10  0 16 17

# After TYPES OF HYPERFACES the type of the hyperface is denoted. As an example, 0 might indicate
# an interior face, 1 a Dirichlet face, and 2 a Neumann face.
# Note that a node can be an interior and a Neumann face, but if it is Dirichlet, all related faces
# have to be of Dirichlet type.

TYPES_OF_HYPERFACES:
0 0 1 1
0 0 1 1
0 0 1 1
0 0 1 1
0 0 1 1
0 0 1 1

# After POINTS_OF_HYPEREDGES the indices of points belonging to a hyperedge have to appear.

POINTS_OF_HYPEREDGES:
 0  1  3  2
 1  4  2  5
 1  6  2  7
 4  8  5 11
 6  8  7 11
 8  9 11 10
