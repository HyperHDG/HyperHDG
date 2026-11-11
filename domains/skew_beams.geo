# Four beams of different lengths and directions in space that meet in one node; the other ends
# are Dirichlet nodes. Two beams start at the common node and two end there, so both orientations
# of a beam against its nodes occur (tests/timoshenko_layouts.cxx).

Space_Dim     = 3;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).

N_Points      = 5;  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = 5;  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  = 4;  # Number of hyperedges that are contained in the hypergraph.

POINTS:
0.0 0.0 0.0
1.0 2.0 0.5
-0.3 0.4 1.1
-2.0 0.1 -0.7
0.2 -0.9 -0.4

HYPERNODES_OF_HYPEREDGES:
0 1
0 2
3 0
4 0

TYPES_OF_HYPERFACES:
0 1
0 1
1 0
1 0

POINTS_OF_HYPEREDGES:
0 1
0 2
3 0
4 0
