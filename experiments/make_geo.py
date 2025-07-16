import matplotlib.pyplot as plt
import numpy as np
import pandas

from datetime import datetime

def make_geo(foldername):
  start_time = datetime.now()
  print("Starting time is", start_time)

  # nodes.csv contains the information about the nodes {n_i}
  # id: The node id <--- removed from the list
  # x,y,z: position in x
  nodes   = pandas.read_csv(foldername + '/nodes.csv')
  nodes   = nodes.to_numpy()[:,1:]
  n_nodes = nodes.shape[0]

  # fibers.csv contains the geometrical information about the fibers
  # id: the fiber id <--- removed from the list
  # node1: node id of node 1
  # node2: node id of node 2
  fibers   = pandas.read_csv(foldername + '/fibers.csv')
  fibers   = fibers.to_numpy()[:,1:]
  n_fibers = fibers.shape[0]

  # connections.csv contains the geometrical information about the connections
  # id: the connection id <--- removed from the list
  # fiber1: fiber id of first edge in connection
  # fiber2: fiber id of second edge in the connection
  # a1: affine constant along fiber1 where the connection is.  
  # a2: affine constant along fiber2 where the connection is.  
  #     ConnectionPos1 = (1-a1)* fiber1.FirstNode() + (a1)* fiber1.SecondNode()
  #     ConnectionPos2 = (1-a2)* fiber2.FirstNode() + (a2)* fiber2.SecondNode()
  connections   = pandas.read_csv(foldername + '/connections.csv')
  connections   = connections.to_numpy()[:,1:]
  n_connections = connections.shape[0]

  # <fiber/connections>Props.csv
  # Id: Id of the edge/connection <--- removed from the list
  # (EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2):   6 structural constants
  # (n_11,n_12,n_13) : normal 1
  # (n_21,n_22,n_23) : normal 2
  fibersProps   = pandas.read_csv(foldername + '/fibersProps.csv')
  fibersProps   = fibersProps.to_numpy()[:,1:]
  n_fibersProps = fibersProps.shape[0]

  connectionsProp   = pandas.read_csv(foldername + '/connectionsProp.csv')
  connectionsProp   = connectionsProp.to_numpy()[:,1:]
  n_connectionsProp = connectionsProp.shape[0]

  end_time = datetime.now()
  print("Reading files ended at", end_time, "after", end_time-start_time)


  vertices, edges, act_fibers = [], [], []
  for con in connections:
    point_a = (1.-con[2]) * nodes[fibers[int(con[0]),0]] + con[2] * nodes[fibers[int(con[0]),1]]
    point_b = (1.-con[3]) * nodes[fibers[int(con[1]),0]] + con[3] * nodes[fibers[int(con[1]),1]]
    point_a, point_b = np.array(point_a), np.array(point_b)

    if len(vertices) == 0:  vertices = np.vstack((point_a, point_b))

    index_a = np.argmin(np.linalg.norm(point_a - vertices, axis=1))
    if np.linalg.norm(point_a - vertices[index_a]) > 1e-10:
      index_a = len(vertices)
      vertices = np.vstack((vertices, point_a))
    index_b = np.argmin(np.linalg.norm(point_b - vertices, axis=1))
    if np.linalg.norm(point_b - vertices[index_b]) > 1e-10:
      index_b = len(vertices)
      vertices = np.vstack((vertices, point_b))

    if index_a != index_b:
      edges.append(np.array([index_a, index_b]))

    act_fibers.append(np.array([con[0], con[2]]))
    act_fibers.append(np.array([con[1], con[3]]))

  end_time = datetime.now()
  print("Preparing data ended at", end_time, "after", end_time-start_time)

  edges_prop = []

  for index in range(n_fibers):
    helper = [ x[1] for x in act_fibers if int(x[0]) == index ]
    if helper == [] or helper is None:  continue

    point_a = nodes[fibers[index,0]]
    point_b = nodes[fibers[index,1]]
    if not any((point_a == x).all() for x in vertices):  vertices = np.vstack((vertices, point_a))
    if not any((point_b == x).all() for x in vertices):  vertices = np.vstack((vertices, point_b))

    helper = list(set(helper))
    helper.sort()
    helper = [0.] + helper + [1.]

    for k in range(len(helper)-1):
      point_ab = (1.-helper[k+0]) * point_a + helper[k+0] * point_b
      point_ba = (1.-helper[k+1]) * point_a + helper[k+1] * point_b

      index_a = np.argmin(np.linalg.norm(point_ab - vertices, axis=1))
      if np.linalg.norm(point_ab - vertices[index_a]) > 1e-10:  print("Error")
      index_b = np.argmin(np.linalg.norm(point_ba - vertices, axis=1))
      if np.linalg.norm(point_ba - vertices[index_b]) > 1e-10:  print("Error")

      if index_a != index_b:
        edges.append(np.array([index_a, index_b]))
        edges_prop.append(fibersProps[index])

  edges_prop = np.vstack((connectionsProp, np.array(edges_prop)))

  min_x, min_y, min_z, max_x, max_y, max_z = 1e10, 1e10, 1e10, -1e10, -1e10, -1e10
  for vertex in vertices:
    min_x, min_y, min_z = min(min_x, vertex[0]), min(min_y, vertex[1]), min(min_z, vertex[2])
    max_x, max_y, max_z = max(max_x, vertex[0]), max(max_y, vertex[1]), max(max_z, vertex[2])

  with open(foldername + '/fiber_network_' + str(len(edges)) + '.geo', 'w') as file:
    file.write("# This file was auto-generated!\n\n")
    file.write("Space_Dim     = 3;  # Dimension of space.\n")
    file.write("HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).\n")
    file.write("N_Points      = " + str(len(vertices)) + ";  # Number of vertices.\n")
    file.write("N_HyperNodes  = " + str(len(vertices)) + ";  # Number of hypernodes.\n")
    file.write("N_HyperEdges  = " + str(len(edges)) + ";  # Number of hyperedges.\n")
    file.write("\nPOINTS:\n")
    for vertex in vertices:
      file.write(str(vertex[0]) + "  " + str(vertex[1]) + "  " + str(vertex[2]) + "\n")
    file.write("\nHYPERNODES_OF_HYPEREDGES:\n")
    for edge in edges:
      file.write(str(int(edge[0])) + "  " + str(int(edge[1])) + "\n")
    file.write("\nTYPES_OF_HYPERFACES:\n")
    for edge in edges:
      left, right = 0, 0
      vertex = vertices[edge[0]]
      if vertex[0] - min_x < 1e-6 * (max_x - min_x) or max_x - vertex[0] < 1e-6 * (max_x - min_x) \
        or vertex[1] - min_y < 1e-6 * (max_y - min_y) or max_y - vertex[1] < 1e-6 * (max_y - min_y):
        left = 1
      vertex = vertices[edge[1]]
      if vertex[0] - min_x < 1e-6 * (max_x - min_x) or max_x - vertex[0] < 1e-6 * (max_x - min_x) \
        or vertex[1] - min_y < 1e-6 * (max_y - min_y) or max_y - vertex[1] < 1e-6 * (max_y - min_y):
        right = 1
      file.write(str(left) + " " + str(right) + "\n")
    file.write("\nPOINTS_OF_HYPEREDGES:\n")
    for edge in edges:
      file.write(str(int(edge[0])) + "  " + str(int(edge[1])) + "\n")
    file.write("\nHYPEREDGE_PROPERTIES: 12\n")
    for edge in edges_prop:
      file.write(str(edge[0]))
      for prop in edge[1:]:
        file.write("  " + str(prop))
      file.write("\n")

  with open(foldername + '/fiber_network_' + str(len(edges)) + '_points.txt', 'w') as file:
    for vertex in vertices:
      file.write(str(vertex[0]) + "  " + str(vertex[1]) + "  " + str(vertex[2]) + "\n")

  ax = plt.figure().add_subplot(111, projection='3d')

  for edge in edges[:n_connections]:
    ax.plot( [vertices[edge[0]][0], vertices[edge[1]][0]], \
             [vertices[edge[0]][1], vertices[edge[1]][1]], \
             [vertices[edge[0]][2], vertices[edge[1]][2]], \
            'rx-' )
  for edge in edges[n_connections:]:
    ax.plot( [vertices[edge[0]][0], vertices[edge[1]][0]], \
             [vertices[edge[0]][1], vertices[edge[1]][1]], \
             [vertices[edge[0]][2], vertices[edge[1]][2]], \
             'bx-' )
  plt.show()

  plt.savefig(foldername + '/graph.png')

  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)


if __name__ == "__main__":
  make_geo('.')
