#!/usr/bin/env python3
"""Minimal transient VTKHDF test file.

3 line segments (4 points), 5 time steps, one scalar field that
varies sinusoidally in time. Mesh is static; PointData is repeated
per step (n_steps * n_points rows), indexed by PointDataOffsets.
"""
import h5py
import numpy as np

n_steps = 5
n_points = 4
n_cells = 3

# Geometry: 4 points along x, 3 line segments connecting them.
points = np.array([[0, 0, 0],
                   [1, 0, 0],
                   [2, 0, 0],
                   [3, 0, 0]], dtype=np.float32)

# VTK_LINE = 3, two points each
connectivity = np.array([0,1, 1,2, 2,3], dtype=np.int64)
offsets      = np.array([0, 2, 4, 6],    dtype=np.int64)  # n_cells + 1
types        = np.array([3, 3, 3],       dtype=np.uint8)

# Time-varying scalar: value at point i, step k = sin(k) + 0.1*i
times = np.linspace(0.0, 1.0, n_steps, dtype=np.float64)
scalar = np.zeros(n_steps * n_points, dtype=np.float32)
for k in range(n_steps):
    for i in range(n_points):
        scalar[k * n_points + i] = np.sin(2 * np.pi * times[k]) + 0.1 * i

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

    # Length n_steps, one entry per step (all equal for static mesh)
    root.create_dataset("NumberOfPoints",          data=np.full(n_steps, n_points, dtype=np.int64))
    root.create_dataset("NumberOfCells",           data=np.full(n_steps, n_cells,  dtype=np.int64))
    root.create_dataset("NumberOfConnectivityIds", data=np.full(n_steps, len(connectivity), dtype=np.int64))

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
