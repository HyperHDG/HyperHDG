#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pandas, logging, argparse, os, sys
import line_profiler

import prin2

def show_network(n_connections, vertices, edges):
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


@line_profiler.profile
def make_geo(input_folder, output_path=".", show=False):
  logger = logging.getLogger("make_geo")

  logger.info("reading files")

  # nodes.csv is the R3 embedding of the ends of the fibers
  #   Id,x,y,z
  # a mapping of the node 'Id' to the position 'x,y,z'
  nodes   = pandas.read_csv(input_folder + "/nodes.csv")
  nodes   = nodes.to_numpy()[:,1:]
  n_nodes = nodes.shape[0]

  # fibers.csv is edge-list of the network, without connections
  #   Id,u,v
  # 'u,v' are the node ids of the ends of the edge 'Id'
  fibers   = pandas.read_csv(input_folder + '/fibers.csv')
  fibers   = fibers.to_numpy()[:,1:]
  n_fibers = fibers.shape[0]

  # connections.csv contains the geometrical information about the connections between fibers
  #   Id,fiber1,fiber2,a1,a2
  # fiber1 is connected to fiber2 with 'Id'.
  # The connection is described by two points, one on each fiber.
  # They can be computed by the mathematical formulae
  #     (1-a1) fiber1.FirstNode() + a1 fiber1.SecondNode(),
  #     (1-a2) fiber2.FirstNode() + a2 fiber2.SecondNode(),
  # where we embed nodes in R3 via the mapping above.
  connections   = pandas.read_csv(input_folder + '/connections.csv')
  connections   = connections.to_numpy()[:,1:]
  n_connections = connections.shape[0]

  # {fiber,connections}Props.csv contain the 6 structural constants
  # associated to each fiber or connection and the 2 normals (principle axis?)
  #   Id, (EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2), (n_11,n_12,n_13), (n_21,n_22,n_23)
  # TODO: why do we need structural
  fibersProps   = pandas.read_csv(input_folder + '/fibersProps.csv')
  fibersProps   = fibersProps.to_numpy()[:,1:]
  n_fibersProps = fibersProps.shape[0]

  connectionsProp   = pandas.read_csv(input_folder + '/connectionsProp.csv')
  connectionsProp   = connectionsProp.to_numpy()[:,1:]
  n_connectionsProp = connectionsProp.shape[0]

  logger.info("preparing data")

  vertices, edges, act_fibers = [], [], []
  for con in connections:
    # con[0..3] = (f1,f2,a1,a2)
    point_a = (1.-con[2]) * nodes[fibers[int(con[0]),0]] + con[2] * nodes[fibers[int(con[0]),1]]
    point_b = (1.-con[3]) * nodes[fibers[int(con[1]),0]] + con[3] * nodes[fibers[int(con[1]),1]]

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

    act_fibers.append([con[0], con[2]])
    act_fibers.append([con[1], con[3]])

  act_fibers = np.array(act_fibers)

  logger.info("creating network")

  if os.path.isdir(output_path):
    output_path += f"/fiber_network_{len(edges)}"
  else:
    # test if we can write to the output path
    with open(output_path + ".geo", "w") as f:
      f.write("test")

  edges_prop = []

  for index in range(n_fibers):
    helper = act_fibers[act_fibers[:,0] == index, 1]
    if helper.size == 0:
      continue

    point_a = nodes[fibers[index,0]]
    point_b = nodes[fibers[index,1]]
    distances_a = np.linalg.norm(point_a - vertices, axis=1)
    if np.min(distances_a) > 1e-15:
      vertices = np.vstack((vertices, point_a))

    distances_b = np.linalg.norm(point_b - vertices, axis=1)
    if np.min(distances_b) > 1e-15:
      vertices = np.vstack((vertices, point_b))


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

  with open(output_path + '.geo', 'w') as file:
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


  with open(output_path + "_points.txt", "w") as file:
    for vertex in vertices:
      file.write(str(vertex[0]) + "  " + str(vertex[1]) + "  " + str(vertex[2]) + "\n")

  if show:
    show_network(n_connections, vertices, edges)


if __name__ == "__main__":
  logging.setLoggerClass(prin2.Logger)
  logger = logging.getLogger("make_geo")

  parser = argparse.ArgumentParser(description="make_geo by Joseph Holten")
  parser.add_argument("folder", help="the folder containing the network data")
  parser.add_argument("-s", "--show", action="store_true",  help="show the created network embedded R3 with plt")
  parser.add_argument("-o", "--output", help="the path where the output is stored", default=".")
  args = parser.parse_args()

  logger.log_args(args)

  make_geo(args.folder, args.output, show=args.show)
