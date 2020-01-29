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

# After "HYPERNODES: the indices of points belonging to a hypernode have to appear. The points need
# to be ordered according to some to be defined scheme and the number of hypernodes has to be equal
# to "Number_of_HyperNodes"! TODO!

HYPERNODES:
0
1
2

# After "HYPEREDGES: the indices of hypernodes belonging to a hyperedge have to appear. These 
# hypernodespoints need to be ordered according to some to be defined scheme and the number of
# hyperedges has to be equal to "Number_of_HyperEdges"! TODO!

HYPEREDGES:
0 1
0 2
1 2
