# --------------------------------------------------------------------------------------------------
# Python functions to write a star_{level}_{dimension}.geo -file
# containing the domain information for HyperHDG
# Ville-Petteri Manninen 28.02.24
# --------------------------------------------------------------------------------------------------
import os
import numpy as np
import math

# --------------------------------------------------------------------------------------------------
# Function for star_{level}_{dimension}.geo -file creation, if one doesn't exist
# --------------------------------------------------------------------------------------------------

def create_star_geo(dimension, level, filename):
  # --- Create a unit cube of dimension (dimension) representing a HyperEdge:

  # Default to 0.0 and 1.0 in 1D
  refPoints = np.array([[0], [1]])
  d = dimension
  newPoints = refPoints

  # Create the other dimensions:
  while d > 1:
    # Repeat the previous dimensions 2 times vertically
    newPoints = np.tile(newPoints, (2, 1))
    # Create the new dimension as (0 ... 0, 1 ... 1)^T
    newDim1 = np.tile(refPoints[0], (int(newPoints.shape[0] / 2), 1))
    newDim2 = np.tile(refPoints[1], (int(newPoints.shape[0] / 2), 1))
    newDim = np.vstack((newDim1, newDim2))

    # Concatenate the new dimension to the old
    newPoints = np.hstack((newPoints, newDim))
    d -= 1
  edgePoints = newPoints

  # --- Create a function to position and scale the HyperEdges (default: origin):

  # Width of the HyperEdges wrt. level
  if (level != 0):
    width = 1 / (2 ** level + 1)
  else:
    width = 1

  # HyperEdge transfomation function to move and scale the HyperEdge
  edgeTransform = lambda position: (edgePoints * width) + position

  # --- Construct a matrix containing all the HyperEdge positions:

  # The range of values a position can take in a dimension
  positions = np.arange(0, 1, width)[:, np.newaxis]
  newPos = positions

  # The median value of the positions
  medVal = positions[math.floor(positions.shape[0]/2)]

  # Again, go through the dimensions and repeat the procedure to include the next one
  d = dimension
  while d > 1:
    # Concatenate the old position with the repeated median in the new dimension:
    topMat = np.hstack((newPos, np.tile(medVal, (newPos.shape[0], 1))))

    # Concatenate the new positions with the repeated median in the old dimensions:
    botMat = np.hstack((np.tile(medVal, (positions.shape[0], newPos.shape[1])), positions))

    # Concatenate, to create the overall positions in the new dimensions
    newPos = np.vstack((topMat, botMat))
    d -= 1

  edgePositions = newPos[np.sort(np.unique(newPos, axis=0, return_index=True)[1]), :]

  # Save the number of HyperEdges:
  N_HyperEdges = edgePositions.shape[0]

  # --- Create the points matrix using the information above:

  # Go through the HyperEdge positionsnp.where(a[:, 0] == 0)
  for i in range(N_HyperEdges):
    # Add the points of the transformed HyperEdge to the overall points:
    if (i == 0):
      points = edgeTransform(edgePositions[i, :])
    else:
      points = np.vstack((points, edgeTransform(edgePositions[i, :])))

  # Take first occurance of each unique point and remove the repeated ones
  points = points[np.sort(np.unique(points, axis=0, return_index=True)[1]), :]

  # --- Points of HyperEdges:

  # Again, go through the HyperEdge positions
  for i in range(N_HyperEdges):
    # Transform the HyperEdge to the position and memorize the points
    edgePointsTf = edgeTransform(edgePositions[i, :])

    # Go through the transformed HyperEdge points
    for j in range(edgePointsTf.shape[0]):
      # Find the index of the point being considered
      rowIndx = np.where((points == edgePointsTf[j, :]).all(axis=1))[0]

      # Concatenate the new index to the other indices to produce a row of indices
      if (j == 0):
        edgeIndx = rowIndx
      else:
        edgeIndx = np.hstack((edgeIndx, rowIndx))

    # Concatenate the indices of the HyperEdge to the overall matrix
    if (i == 0):
      pointsOfHyEdges = edgeIndx
    else:
      pointsOfHyEdges = np.vstack((pointsOfHyEdges, edgeIndx))

  # --- Hypernodes of Hyperedges:

  # Construct a matrix to reference the points of the HyperEdge to Nodes:

  # Go through the dimensions
  for i in range(dimension):
    # And the point choices
    for j in refPoints:
      # Find the point indices on the referece edge for the node
      newNode = np.where(edgePoints[:, i] == j)
      try:
        refNodePoints = np.vstack((refNodePoints, newNode))
      except NameError:
        refNodePoints = newNode

  # Construct the a matrix to contain the indices of the Points of the HyperNodes:

  # Go through the HyperEdges
  for i in range(N_HyperEdges):
    # Go through the reference point indices for the HyperNodes
    for j in range(refNodePoints.shape[0]):
      # Using the reference indices, find the correct point indices
      if (level == 0):
        # Level = 0, has lower dimension of matrices:
        newNodePoints = pointsOfHyEdges[refNodePoints[j, :]]
      else:
        newNodePoints = pointsOfHyEdges[i, refNodePoints[j, :]]

      # Add the new row of Points to the Points of HyperNodes matrix
      try:
        pointsOfNodes = np.vstack((pointsOfNodes, newNodePoints))
      except NameError:
        pointsOfNodes = newNodePoints

  # Sort the rows and take unique instances for later:
  sPointsOfNodes = np.sort(pointsOfNodes, axis=1)
  pointsOfNodes = sPointsOfNodes[np.sort(np.unique(sPointsOfNodes, axis=0, return_index=True)[1]), :]

  # Construct the Nodes of HyperEdges matrix

  # Go through the HyperEdges
  for i in range(N_HyperEdges):
    # ... and the reference point indices (again)
    for j in range(refNodePoints.shape[0]):
      # Using the reference indices, find the correct point indices
      if (level == 0):
        nodePoints = pointsOfHyEdges[refNodePoints[j, :]]
      else:
        nodePoints = pointsOfHyEdges[i, refNodePoints[j, :]]

      # Find the correct HyperNode using the point information
      newNode = np.where((pointsOfNodes == np.sort(nodePoints)).all(axis=1))[0]

      # ... and add it into the HyperEdge's nodes:
      if (j == 0):
        newEdge = newNode
      else:
        newEdge = np.hstack((newEdge, newNode))

    # Add the HyperEdge to the overall matrix:
    if (i == 0):
      nodesOfHyEdges = newEdge
    else:
      nodesOfHyEdges = np.vstack((nodesOfHyEdges, newEdge))

  # --- Finally, construct the Types of HyperFaces matrix

  # Find the HyperEdges on the boundary
  boundaryEdges = np.argwhere((edgePositions == np.min(edgePositions)) | (edgePositions == np.max(edgePositions)))

  # Build the matrix itself:
  if (level == 0):
    # Initialize the types as matrix of ones:
    typesOfHyFaces = np.ones(nodesOfHyEdges.shape, dtype=int)
  else:
    # Initialize the types as matrix of zeros:
    typesOfHyFaces = np.zeros(nodesOfHyEdges.shape, dtype=int)
    # Set the correct node as the boundary node for the edges on the boundary:
    for indxs in boundaryEdges:
      if (edgePositions[indxs[0], indxs[1]] == np.min(edgePositions)):
        typesOfHyFaces[indxs[0], 2 * indxs[1]] = 1
      else:
        typesOfHyFaces[indxs[0], 2 * indxs[1] + 1] = 1

  # --- Write the file itself
  N_points = points.shape[0]
  N_HyperNodes = np.max(nodesOfHyEdges + 1)

  f = open(os.path.dirname(os.path.abspath(__file__)) + \
           "/../domains/" + filename, 'w') # Create and open the file for writing

  # Write the information into the .geo -file
  f.write(f'''# File containing geometric and topological information for of hypergraph being an "star shape" of level {level} and dimension {dimension}.
        
Space_Dim     = {dimension};  # Dimension of space indicating amount of coordinates of points.
HyperEdge_Dim = {dimension};  # Dimension of hyperedge (must be uniform).

N_Points      =  {N_points};  # Number of vertices that are contained in the hypergraph.
N_HyperNodes  = {N_HyperNodes};  # Number of hypernodes that are contained in the hypergraph.
N_HyperEdges  =  {N_HyperEdges};  # Number of hyperedges that are contained in the hypergraph.

POINTS:\n''')

  if (level == 0):
    for i in range(N_points):
      for j in range(dimension):
        f.write(f'{points[i, j]} ')
      f.write('\n')

    f.write('\nHYPERNODES_OF_HYPEREDGES:\n')

    for i in range(N_HyperEdges):
      for j in range(nodesOfHyEdges.shape[0]):
        f.write(f'{nodesOfHyEdges[j]} ')
      f.write('\n')

    f.write('\nTYPES_OF_HYPERFACES:\n')

    for i in range(N_HyperEdges):
      for j in range(typesOfHyFaces.shape[0]):
        f.write(f'{typesOfHyFaces[j]} ')
      f.write('\n')

    f.write('\nPOINTS_OF_HYPEREDGES:\n')

    for i in range(N_HyperEdges):
      for j in range(pointsOfHyEdges.shape[0]):
        f.write(f'{pointsOfHyEdges[j]} ')
      f.write('\n')

  else:
    for i in range(N_points):
      for j in range(dimension):
        f.write(f'{points[i, j]} ')
      f.write('\n')

    f.write('\nHYPERNODES_OF_HYPEREDGES:\n')

    for i in range(N_HyperEdges):
      for j in range(nodesOfHyEdges.shape[1]):
        f.write(f'{nodesOfHyEdges[i, j]} ')
      f.write('\n')

    f.write('\nTYPES_OF_HYPERFACES:\n')

    for i in range(N_HyperEdges):
      for j in range(typesOfHyFaces.shape[1]):
        f.write(f'{typesOfHyFaces[i, j]} ')
      f.write('\n')

    f.write('\nPOINTS_OF_HYPEREDGES:\n')

    for i in range(N_HyperEdges):
      for j in range(pointsOfHyEdges.shape[1]):
        f.write(f'{pointsOfHyEdges[i, j]} ')
      f.write('\n')
    return

# --------------------------------------------------------------------------------------------------
# Function to get the path of the defined star_{level}.geo -file
# --------------------------------------------------------------------------------------------------

def get_star_domain(dimension, level):
  filename = 'star_' + str(dimension) + '_' + str(level) + '.geo'
  if (os.path.exists(os.path.dirname(os.path.abspath(__file__)) + \
      "/../domains/" + filename)):
    return filename

  create_star_geo(dimension, level, filename)
  return filename
