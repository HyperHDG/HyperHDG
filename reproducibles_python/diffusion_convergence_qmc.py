import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from joblib import Parallel, delayed

import os, sys

try:
	import HyperHDG
except(ImportError, ModuleNotFoundError) as error:
	sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
	import HyperHDG

pol_g = 1

n_shifts = 16
gen_vec = HyperHDG.qmc_methods.generating_vector("lattice-32001-1024-1048576.3600", 100)


zfgen = np.random.default_rng(27)

const = HyperHDG.config()
const.global_loop = "Elliptic"
const.topology = "Cubic<2, 2>"
const.geometry = "UnitCube<2, 2, double>"
const.node_descriptor = "Cubic<2, 2>"
const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionAffine, double, std::vector<double> >"
const.param_t = "std::vector<double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", "double", "double", "vector[double]"]
const.include_files   = ["reproducibles_python/parameters/diffusion_qmc.hxx"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP([100, 100])

system_size = HDG_wrapper.size_of_system()
mean = np.zeros(system_size)
mean_sa = np.zeros(system_size)

with open("output/kvgr.txt", "w") as f:
	f.write("n_qmc_points\tsq_err\n")

def solve(k, gen_vec, shift, n_qmc_points, HDG_wrapper):
	print("Punkt " + str(k+1) + " von " + str(n_qmc_points))
	system_size = HDG_wrapper.size_of_system()
	arr = HyperHDG.qmc_methods.get_quadrature_point(gen_vec, shift, k+1, n_qmc_points) - 0.5
	vectorRHS = np.multiply( HDG_wrapper.residual_flux(np.zeros(system_size), arr), -1. )
	col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(arr)
	A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
	[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-11)
	return HDG_wrapper.errors(vectorSolution)

for i in range(2, 4):
	n_qmc_points = 2**i
	print("Quadraturpunkte: " + str(n_qmc_points))
	mean_sa = 0
	for m in range(n_shifts):
		mean_sa = 0
		print("Shift " + str(m+1) + " von " + str(n_shifts))
		shift = zfgen.random(100)
		errs = Parallel(n_jobs=-2, prefer="threads")( delayed(solve) (k, gen_vec, shift, n_qmc_points, HDG_wrapper) for k in range(n_qmc_points))
		mean = np.mean(np.array(errs))
		mean_sa += mean / n_shifts
	with open("output/kvgr.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(np.linalg.norm(mean_sa - mean)) + "\n")

