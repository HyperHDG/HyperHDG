#!/usr/bin/env python3
"""LDG-H diffusion on the unit square from a module that CMake built in advance.

    cmake -S . -B build && cmake --build build --target diffusion_elliptic_nb
    PYTHONPATH=build/examples python3 examples/diffusion_elliptic_nb_precompiled.py

The module is examples/diffusion_elliptic_nb_module.cxx: one instantiation (dimension 2, degree 1)
of GlobalLoop::Generic with LocalSolver::Diffusion. examples/diffusion_elliptic_nb_jit.py does the
same with any dimension and degree, compiled on first use.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import diffusion_elliptic_nb

n = 8  # elements per direction
hdg = diffusion_elliptic_nb.Diffusion([n, n])
hdg.constrain_dirichlet_nodes()
hdg.prepare()

# The system: the residual at zero is minus the right-hand side, assemble writes the matrix
# (a unit diagonal on the Dirichlet dofs) into arrays sized by n_triplets.
size = hdg.size_of_system()
rhs = np.empty(size)
hdg.residual(np.zeros(size), rhs)
rhs *= -1.
m = hdg.n_triplets(True)
rows, cols, vals = np.empty(m, np.uint32), np.empty(m, np.uint32), np.empty(m)
hdg.assemble(rows, cols, vals, 0., True)
matrix = sp.coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsc()

solution = spla.spsolve(matrix, rhs)
print(f"unknowns: {size}")
print(f"error: {hdg.errors(solution)[0]:.6e}")
