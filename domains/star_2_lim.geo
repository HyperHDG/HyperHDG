# File containing geometric and topological information for of hypergraph being a limit +-shape.
        
Space_Dim     = 2;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).

N_Points      = 5;  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = 5;  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  = 4;  # Number of hyperedges that are contained in the hypergraph.

POINTS:
0.0 0.5
0.5 0.5
1.0 0.5
0.5 0.0
0.5 1.0

HYPERNODES_OF_HYPEREDGES:
0 1
1 2
3 1
1 4

TYPES_OF_HYPERFACES:
1 0
0 1
1 0
0 1

POINTS_OF_HYPEREDGES:
0 1
1 2
3 1
1 4
