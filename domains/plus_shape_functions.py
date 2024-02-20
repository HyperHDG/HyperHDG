# --------------------------------------------------------------------------------------------------
# Python functions to write/read a plus_shape_{level}.geo -file for domain information for HyperHDG
# Ville-Petteri Manninen 19.02.24
# --------------------------------------------------------------------------------------------------

import os, sys
import numpy as np
import math

# --------------------------------------------------------------------------------------------------
# Function for plus_domain_{level}.geo -file creation, if one doesn't exist
# --------------------------------------------------------------------------------------------------
def create_plus_domain_geo(level, filename):
    f = open(os.path.dirname(os.path.abspath(__file__)) + \
             "/../domains/" + filename, 'w')  # Create and open the file for writing
    N_edges = 2 ** (level + 1) + 1  # Number of hyperedges wrt. level
    plus_width = 1 / (math.ceil(N_edges / 2))  # Width of arms of the +-sign
    points = np.zeros((1, 2))  # Initialize the points matrix

    # -- Initialize the vectors for the point coordinates
    # Vertical arms:
    x1 = np.array([0.5 - plus_width / 2, 0.5 + plus_width / 2])
    y1 = np.arange(0, 1 + plus_width / 2, plus_width)

    # Horizontal arms:
    # x and y interchanged, with 4 center points skipped

    # Generate the POINTS -variable as a matrix
    for y in y1:
        for x in x1:
            points = np.vstack([points, [x, y]])

    for x in y1[np.arange(0, math.floor(y1.size / 2) - 1, 1)]:
        for y in x1:
            points = np.vstack([points, [x, y]])

    for x in y1[np.arange(math.ceil(y1.size / 2) + 1, y1.size, 1)]:
        for y in x1:
            points = np.vstack([points, [x, y]])

    points = np.delete(points, 0, 0)
    N_points = points.shape[0]

    # Generate the HYPERNODES_OF_HYPEREDGES -variable as a matrix
    nodes_of_edges = np.zeros((N_edges, 4), dtype=int)
    middle_edge = math.floor(N_edges / 4)

    for i in range(N_edges):
        # First edge has the first 4 nodes:
        if (i == 0):
            nodes_of_edges[i, :] = [i, i + 1, i + 2, i + 3]
            continue

        # Vertical arm from bottom to top:
        if (i < N_edges / 2):
            k = nodes_of_edges[i - 1, 3]
            nodes_of_edges[i, :] = [k + 1, k + 2, k, k + 3]

            # Memorize the middle edge's left and right node:
            if (i == middle_edge):
                left_node = k + 1
                right_node = k + 2
            continue

        # This pattern works for all levels to recognize the edges adjacent to the middle ones:
        if (i == middle_edge * 3):
            # Level 1 has different reference (top most edge):
            if (level == 1):
                k = nodes_of_edges[i - 1, 3]
                nodes_of_edges[i, :] = [k + 1, left_node, k + 2, k + 3]
                continue

            # While the others reference the previous edge to the left:
            k = nodes_of_edges[i - 1, 1]
            nodes_of_edges[i, :] = [k, left_node, k + 3, k + 4]
            continue

        if (i == middle_edge * 3 + 1):
            k = nodes_of_edges[i - 1, 3]
            nodes_of_edges[i, :] = [right_node, k + 1, k + 2, k + 3]
            continue

        # The rest of the edges:
        if (i == math.ceil(N_edges / 2)):
            k = nodes_of_edges[i - 1, 3]
            nodes_of_edges[i, :] = [k + 1, k + 2, k + 3, k + 4]
            continue

        k = nodes_of_edges[i - 1, 1]
        nodes_of_edges[i, :] = [k, k + 3, k + 4, k + 5]

    N_nodes = np.max(nodes_of_edges + 1)

    # Generate the TYPES_OF_HYPERFACES -variable as a matrix
    temp = np.zeros((N_edges, 4), dtype=int)
    M1 = np.where(nodes_of_edges != 2, temp, 1)
    M2 = np.where(nodes_of_edges != N_nodes - 3, temp, 1)
    M3 = np.where(nodes_of_edges != (N_nodes + 2) / 2, temp, 1)
    M4 = np.where(nodes_of_edges != (N_nodes + 2) / 2 + 1, temp, 1)
    types_of_faces = M1 + M2 + M3 + M4

    # Generate the POINTS_OF_HYPEREDGES: -variable as a matrix
    points_of_edges = np.zeros((N_edges, 4), dtype=int)

    for i in range(N_edges):
        # Vertical arm:
        if (i < N_edges / 2):
            k = 2 * i
            points_of_edges[i, :] = [k, k + 1, k + 2, k + 3]
            continue

        # Edge left of the middle one:
        if (i == 3 * middle_edge):
            p = points_of_edges[middle_edge, 0]
            if (level == 1):
                k = points_of_edges[i - 1, 3]
            else:
                k = points_of_edges[i - 1, 2]

            points_of_edges[i, :] = [k + 1, p, k + 2, p + 2]
            continue

        # Edge right of the middle one:
        if (i == 3 * middle_edge + 1):
            k = points_of_edges[i - 1, 2]
            p = points_of_edges[middle_edge, 1]
            points_of_edges[i, :] = [p, k + 1, p + 2, k + 2]
            continue

        # The rest of the horizontal edges:
        if (i == math.ceil(N_edges / 2)):
            k = points_of_edges[i - 1, 3]
            points_of_edges[i, :] = [k + 1, k + 3, k + 2, k + 4]
            continue

        if (i < 3 * middle_edge):
            k = points_of_edges[i - 1, 0]
            points_of_edges[i, :] = [k, k + 3, k + 2, k + 4]

        k = points_of_edges[i - 1, 1]
        points_of_edges[i, :] = [k, k + 2, k + 1, k + 3]

    # Write the information into the .geo -file
    f.write(
        f'''# File containing geometric and topological information for of hypergraph being a +-shape of level {level}.
        
Space_Dim     = 2;  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = 2;  # Dimension of hyperedge (must be uniform).

N_Points      =  {N_points};  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = {N_nodes};  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  =  {N_edges};  # Number of hyperedges that are contained in the hypergraph.

POINTS:\n''')
    for i in range(N_points):
        f.write(f'{points[i, 0]}\t {points[i, 1]}\n')

    f.write('\nHYPERNODES_OF_HYPEREDGES:\n')

    for i in range(N_edges):
        for j in range(4):
            f.write(f'{nodes_of_edges[i, j]}\t')

        f.write('\n')

    f.write('\nTYPES_OF_HYPERFACES:\n')

    for i in range(N_edges):
        for j in range(4):
            f.write(f'{types_of_faces[i, j]}\t')

        f.write('\n')

    f.write('\nPOINTS_OF_HYPEREDGES:\n')

    for i in range(N_edges):
        for j in range(4):
            f.write(f'{points_of_edges[i, j]}\t')

        f.write('\n')

    return


# --------------------------------------------------------------------------------------------------
# Function to fetch the path of the defined plus_domain_{level}.geo -file
# --------------------------------------------------------------------------------------------------
def get_plus_domain(level):
    filename = 'plus_domain_' + str(level) + '.geo'

    if (os.path.exists(os.path.dirname(os.path.abspath(__file__)) + \
                       "/../domains/" + filename)):
        return filename

    create_plus_domain_geo(level, filename)
    return filename