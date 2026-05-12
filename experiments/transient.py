#!/usr/bin/env python3
"""Minimal transient VTKHDF test file.

3 line segments (4 points), 5 time steps, one scalar field that
varies sinusoidally in time. Mesh is static; PointData is repeated
per step (n_steps * n_points rows), indexed by PointDataOffsets.
"""
import h5py
import numpy as np

n_steps = 20
n_points = 20
n_cells = n_points - 1

# Geometry: linspace between 0 and 1 along x.
points = np.zeros((n_points, 3), dtype=np.float32)
points[:, 0] = np.linspace(0, 1, n_points, dtype=np.float32)

# VTK_LINE = 3, two points each, chained
connectivity = np.empty(2 * n_cells, dtype=np.int64)
connectivity[0::2] = np.arange(n_cells)
connectivity[1::2] = np.arange(1, n_cells + 1)
offsets = np.arange(n_cells + 1, dtype=np.int64) * 2
types   = np.full(n_cells, 3, dtype=np.uint8)

# Field: sin(2*pi*x) * cos(t), t in [0, 2*pi]
times = np.linspace(0.0, 2*np.pi, n_steps, dtype=np.float64)
scalar = np.zeros(n_steps * n_points, dtype=np.float32)
for k in range(n_steps):
    for i in range(n_points):
        scalar[k * n_points + i] = np.sin(2 * np.pi * points[i, 0]) * np.cos(times[k])

# PointDataOffsets: row in PointData/scalar where each step starts
pd_offsets = np.arange(n_steps, dtype=np.int64) * n_points

with h5py.File("transient.vtkhdf", "w") as f:
    root = f.create_group("VTKHDF")
    root.attrs.create("Version", [2, 0], dtype=np.int64)
    root.attrs.create("Type", np.bytes_("UnstructuredGrid"))

    # Geometry (static, single-part: NumberOf* are length-1 arrays)
    root.create_dataset("Points", data=points)
    root.create_dataset("Connectivity", data=connectivity)
    root.create_dataset("Offsets", data=offsets)
    root.create_dataset("Types", data=types)

    # NumberOf* as length-1 arrays (single part, static mesh)
    root.create_dataset("NumberOfPoints",          data=np.array([n_points], dtype=np.int64))
    root.create_dataset("NumberOfCells",           data=np.array([n_cells],  dtype=np.int64))
    root.create_dataset("NumberOfConnectivityIds", data=np.array([len(connectivity)], dtype=np.int64))

    # PointData
    pd = root.create_group("PointData")
    pd.create_dataset("scalar", data=scalar)

    # Steps
    steps = root.create_group("Steps")
    steps.attrs.create("NSteps", n_steps, dtype=np.int64)
    steps.create_dataset("Values", data=times)
    steps.create_dataset("PartOffsets", data=np.zeros(n_steps, dtype=np.int64))
    steps.create_dataset("PointOffsets",         data=np.zeros(n_steps, dtype=np.int64))
    steps.create_dataset("CellOffsets",          data=np.zeros(n_steps, dtype=np.int64))
    steps.create_dataset("ConnectivityIdOffsets", data=np.zeros(n_steps, dtype=np.int64))
    steps.create_dataset("NumberOfParts", data=np.ones(n_steps, dtype=np.int64))

    pdo = steps.create_group("PointDataOffsets")
    pdo.create_dataset("scalar", data=pd_offsets)

print("wrote transient.vtkhdf")
