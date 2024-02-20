# --------------------------------------------------------------------------------------------------
# Python functions to write a Lshape_{level}.geo -file to read domain information from for HyperHDG
# Ville-Petteri Manninen 14.02.24
# --------------------------------------------------------------------------------------------------


import os, sys
import numpy as np

# --------------------------------------------------------------------------------------------------
# Function for Lshape_{level}.geo -file creation, if one doesn't exist
# --------------------------------------------------------------------------------------------------

def create_l_domain_geo(level, filename):
  f = open(os.path.dirname(os.path.abspath(__file__)) + \
    "/../domains/" + filename, 'w')  # Create and open the file for writing
  L_width = 1 / 2 ** level  # Width of the L-shape wrt. level
  points = np.zeros((1, 2))  # Initialize the points matrix
  N_edges = 2 ** (level + 1) - 1  # Number of hyperedges wrt. level

  # Initialize the vectors for the point coordinates
  y1 = np.arange(0, 1 + L_width, L_width)
  x1 = np.array([0, L_width])
  y2 = np.array([1 - L_width, 1])
  x2 = np.arange(2 * L_width, 1 + L_width, L_width)

  # Generate the POINTS -variable as a matrix

  for y in y1:
    for x in x1:
      points = np.vstack([points, [x, y]])

  for x in x2:
    for y in y2:
      points = np.vstack([points, [x, y]])

  points = np.delete(points, 1, 0)
  N_points = points.shape[0]

  # Generate the HYPERNODES_OF_HYPEREDGES -variable as a matrix

  nodes_of_edges = np.zeros((N_edges, 4), dtype=int)

  for i in range(N_edges):
    if (i == 0):
      nodes_of_edges[i, :] = [i, i + 1, i + 2, i + 3]
      continue

    if (i < N_edges / 2):
      k = nodes_of_edges[i - 1, 3]
      nodes_of_edges[i, :] = [k + 1, k + 2, k, k + 3]
      continue

    if (i < N_edges / 2 + 1):
      k = nodes_of_edges[i - 1, 1]
      nodes_of_edges[i, :] = [k, k + 2, k + 3, k + 4]
      continue

    k = nodes_of_edges[i - 1, 1]
    nodes_of_edges[i, :] = [k, k + 3, k + 4, k + 5]

  N_nodes = np.max(nodes_of_edges + 1)

  # Generate the TYPES_OF_HYPERFACES -variable as a matrix
  temp = np.zeros((N_edges, 4), dtype=int)
  M1 = np.where(nodes_of_edges != 2, temp, 1)
  M2 = np.where(nodes_of_edges != N_nodes - 3, temp, 1)
  types_of_faces = M1 + M2

  # Generate the POINTS_OF_HYPEREDGES: -variable as a matrix

  points_of_edges = np.zeros((N_edges, 4), dtype=int)
  for i in range(N_edges):
    if (i < N_edges / 2):
      k = 2 * i
      points_of_edges[i, :] = [k, k + 1, k + 2, k + 3]
      continue

    if (i < N_edges / 2 + 1):
      k = points_of_edges[i - 1, 1]
      points_of_edges[i, :] = [k, k + 3, k + 2, k + 4]
      continue

    k = points_of_edges[i - 1, 1]
    points_of_edges[i, :] = [k, k + 2, k + 1, k + 3]

  # Write the information into the .geo -file
  f.write(f'''# File containing geometric and topological information for of hypergraph being an L-shape of level {level}.

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
      f.write(f'{nodes_of_edges[i, j]} ')

    f.write('\n')

  f.write('\nTYPES_OF_HYPERFACES:\n')

  for i in range(N_edges):
    for j in range(4):
      f.write(f'{types_of_faces[i, j]} ')

    f.write('\n')

  f.write('\nPOINTS_OF_HYPEREDGES:\n')

  for i in range(N_edges):
    for j in range(4):
      f.write(f'{points_of_edges[i, j]} ')

    f.write('\n')

  return


# --------------------------------------------------------------------------------------------------

# Function to get the path of the defined Lshape_{level}.geo -file

# --------------------------------------------------------------------------------------------------

def get_l_domain(level):
  filename = 'L_domain_' + str(level) + '.geo'

  if (os.path.exists(os.path.dirname(os.path.abspath(__file__)) + \
      "/../domains/" + filename)):
    return filename

  create_l_domain_geo(level, filename)
  return filename
