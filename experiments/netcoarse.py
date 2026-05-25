#!/usr/bin/env python
"""
netcoarse.py — convert a PETSc matrix to a transient VTKHDF file.

Each column of the input matrix becomes one time step of a scalar PointData
field 'values' on a network of line cells (VTK_LINE) defined by the domain
file. Useful for visualizing basis functions, eigenmodes, or column-indexed
fields over a beam/fiber network in ParaView.

Open in ParaView GUI or visualize with netvis.py via

experiments/netvis.py <VTKHDF> -r <RADIUS> --view iso --warp-scale <SCALE> \
  --arrows 0 --warp-by values:0:kHat

Inputs
------
domain : HDF5 file with
    /domain/points  (Np, 3) float  -- node coordinates
    /domain/edges   (Ne, 2) int    -- point indices, one edge per cell
matrix : PETSc binary Mat of shape (Np, Ncols). Rows index domain points,
    columns index the fields/modes to visualize as time steps.

Output
------
VTKHDF v2.0 UnstructuredGrid. Geometry is shared across all steps; only the
PointData window into "values" advances Np per step.

Requirements
------------
PETSC_DIR must be set; PetscBinaryIO is loaded from $PETSC_DIR/lib/petsc/bin.

Caveats
-------
The matrix is densified in memory (A.todense). Will OOM for very large
networks or column counts; rewrite with a streamed column loop if needed.

Author
------
Joseph Holten, KIT, 2026.
"""

import h5py, numpy as np, sys, os
sys.path.append(os.path.join(os.environ['PETSC_DIR'], "lib/petsc/bin"))
import PetscBinaryIO
import scipy.sparse as sp
import argparse

def main(domain_path, mat_path, out_path):
  # --- read PETSc matrix: rows = nodes, cols = basis functions
  io = PetscBinaryIO.PetscBinaryIO()
  fh = open(mat_path)
  otype = io.readObjectType(fh)
  assert otype == "Mat"
  M = io.readMat(fh, mattype="scipy.sparse")
  A = np.asarray(M.todense(), dtype=np.float32)
  n_pts, n_cols = A.shape

  with h5py.File(domain_path, 'r') as d:
    points = d['/domain/points'][:].astype(np.float32)             # (Np, 3)
    edges  = d['/domain/edges'][:].astype(np.int64)                # (Ne, 2)

  assert points.shape[0] == n_pts, \
    f"matrix rows ({n_pts}) must equal domain points ({points.shape[0]})"

  n_cells = edges.shape[0]
  connectivity = edges.reshape(-1).astype(np.int64)
  offsets = (np.arange(n_cells + 1) * 2).astype(np.int64)
  types = np.full(n_cells, 3, dtype=np.uint8) # VTK_LINE = 3

  vals = A.T.reshape(-1).astype(np.float32) # (n_cols * n_pts,)

  # --- write VTKHDF
  with h5py.File(out_path, 'w') as f:
    g = f.create_group('VTKHDF')
    dt = h5py.string_dtype(encoding='ascii', length=16)
    g.attrs.create('Type', 'UnstructuredGrid', dtype=dt)
    g.attrs['Version'] = np.array([2, 0], dtype=np.int64)            # check

    g.create_dataset('Points', data=points, maxshape=(None, 3))
    g.create_dataset('Connectivity', data=connectivity, maxshape=(None,))
    g.create_dataset('Offsets', data=offsets, maxshape=(None,))
    g.create_dataset('Types', data=types, maxshape=(None,))
    g.create_dataset('NumberOfPoints', data=np.array([n_pts], dtype=np.int64),
                     maxshape=(None,))
    g.create_dataset('NumberOfCells', data=np.array([n_cells], dtype=np.int64),
                     maxshape=(None,))
    g.create_dataset('NumberOfConnectivityIds',
                     data=np.array([connectivity.size], dtype=np.int64),
                     maxshape=(None,))

    pd = g.create_group('PointData')
    pd.create_dataset('values', data=vals, maxshape=(None,))

    s = g.create_group('Steps')
    s.attrs['NSteps'] = np.int64(n_cols)
    zeros = np.zeros(n_cols, dtype=np.int64)
    s.create_dataset('PointOffsets', data=zeros, maxshape=(None,))
    s.create_dataset('CellOffsets', data=zeros, maxshape=(None,))
    s.create_dataset('ConnectivityIdOffsets', data=zeros, maxshape=(None,))
    s.create_dataset('PartOffsets', data=zeros, maxshape=(None,))
    s.create_dataset('NumberOfParts', data=np.ones(n_cols, dtype=np.int64), maxshape=(None,))
    s.create_dataset('Values', data=np.arange(n_cols, dtype=np.float64), maxshape=(None,))
    pdo = s.create_group('PointDataOffsets')
    pdo.create_dataset('values', data=(np.arange(n_cols) * n_pts).astype(np.int64),
                       maxshape=(None,))

if __name__ == "__main__":
  p = argparse.ArgumentParser(
      description="Convert a PETSc Mat (Np x Ncols) into a transient VTKHDF "
                  "UnstructuredGrid over a line network. Each matrix column "
                  "becomes one time step of PointData/values.")
  p.add_argument("domain", help="HDF5 domain file with /domain/points and /domain/edges")
  p.add_argument("matrix", help="PETSc binary Mat; rows must match number of domain points")
  p.add_argument("output", help="Output .vtkhdf path (HDF5, overwritten)")
  a = p.parse_args()
  main(a.domain, a.matrix, a.output)
